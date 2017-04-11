#!/nfs/goldstein/software/python2.7/bin/python2.7.7
"""
Create a configuration file used by the Dragen system to processed IGM samples
based on their priority level in the Dragen pipeline queue
"""

import argparse
import MySQLdb
import os
import sys
import subprocess
from create_align_config import create_align_config
from ConfigParser import SafeConfigParser
from datetime import datetime
from dragen_sample import dragen_sample
from glob import glob

def main(samples, debug, dontexecute, database):
    CNF = "/nfs/goldstein/software/dragen/dragen.cnf" # defaults file for pipeline
    config_parser = SafeConfigParser()
    config_parser.read(CNF)
    parameters = {}
    for parameter in ("DB_CONFIG_FILE","DRAGEN_CONFIG"):
        parameters[parameter] = config_parser.get("dragen_pipe", parameter.upper())

    if debug:
        print "Parameters used: {0}".format(parameters)

    #Automated run
    if samples == True:
        db = get_cursor(parameters,database)
        curs = db.cursor()
        info = get_next_sample(curs,debug)
        if info is None:
            print "No samples were found"
            sys.exit()

        while info is not None:
            dragen_id = info[4]
            db = get_cursor(parameters,database)
            curs = db.cursor()
            sample = dragen_sample(info[0],info[1],info[2],info[3],curs)
            single_sample_setup(curs,sample,parameters,debug)
            db.close()

            run_sample(sample,dragen_id,dontexecute,parameters,database,debug)

            db = get_cursor(parameters,database)
            curs = db.cursor()
            info = get_next_sample(curs,debug)
            db.close()

def get_cursor(parameters,database):
    db = MySQLdb.connect(db=database, read_default_group="clientsequencedb",
            read_default_file=parameters["DB_CONFIG_FILE"])
    return db

def run_sample(sample,dragen_id,dontexecute,parameters,database,debug):
    output_dir = sample.metadata['output_dir']
    pseudo_prepid = sample.metadata['pseudo_prepid']
    sample_bam = "{output_dir}/{sample_name}.bam".format(
            output_dir=output_dir,
            sample_name=sample.metadata['sample_name'])
    db = get_cursor(parameters,database)
    curs = db.cursor()
    submitTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    if os.path.isfile(sample_bam) == False:
        cmd = ['dragen', '-f', '-v', '-c', sample.metadata['conf_file'], '--watchdog-active-timeout', '30000']
        if debug:
            print ' '.join(cmd)

        dragen_stderr = open(sample.metadata['dragen_stderr'],'a')
        if dontexecute == False:
            if dragen_id != 0: # All manually run samples will have a dragen_id of 0
                update_queue(curs,dragen_id,debug)

            curs.close()

            with open(sample.metadata['dragen_stdout'],'a') as dragen_stdout:

                process = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=dragen_stderr)
                for line in iter(process.stdout.readline, ''):
                    if debug:
                        sys.stdout.write(line)
                    dragen_stdout.write(line)

                process.communicate()
                error_code = process.wait()
                subprocess.call(['chmod','-R','775','{}'.format(output_dir)])
            if debug:
               print "Dragen error code: {0}".format(error_code)
            if error_code == 0:
                db = get_cursor(parameters,database)
                curs = db.cursor()
                finishTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

                insert_query = ("INSERT INTO gatk_queue "
                    "SELECT * FROM tmp_dragen WHERE dragen_id={0}"
                    ).format(dragen_id)
                rm_query = "DELETE FROM {0} WHERE dragen_id={1}".format("tmp_dragen",dragen_id)

                if debug:
                    print insert_query
                    print rm_query
                update_pipeline_step_id(curs,sample,submitTime,finishTime,debug)
                update_dragen_metadata(curs,sample)
                curs.execute(insert_query)
                curs.execute(rm_query)
                db.commit()
                db.close()

            else:
                print "Sample was not deleted from tmp_dragen due to error code: {0}".format(error_code)

            dragen_stdout.close()
            dragen_stderr.close()
        else:
            print "Sample {sample_name} config complete.".format(**sample.metadata)
    else:
        print "Sample {sample_name} bam file already exists!".format(sample_name=sample.metadata['sample_name'])
        update_queue(curs,dragen_id,debug)
        db.close()

def update_dragen_metadata(curs,sample):
    sample_name = sample.metadata['sample_name']
    pseudo_prepid = sample.metadata['pseudo_prepid']
    sample_type = sample.metadata['sample_type']
    capture_kit = sample.metadata['capture_kit']
    prepid = sample.metadata['prepid'][0]
    priority = sample.metadata['priority']

    query = ("SELECT * from dragen_sample_metadata "
            "WHERE pseudo_prepid = {} "
            ).format(pseudo_prepid)
    curs.execute(query)
    dbInfo = curs.fetchone()
    if dbInfo:
        pass
    else:
        insertQuery = ("INSERT INTO dragen_sample_metadata "
                "(sample_name,pseudo_prepid,sample_type,capture_kit,prepid,priority) "
                "VALUES ('{}',{},'{}','{}',{},{} ) "
                ).format(sample_name,pseudo_prepid,sample_type,capture_kit,prepid,priority)
        print insertQuery
        curs.execute(insertQuery)

def update_pipeline_step_id(curs,sample,submitTime,finishTime,debug):
    pseudo_prepid = sample.metadata['pseudo_prepid']
    timesRunQuery = ("SELECT times_ran FROM dragen_pipeline_step "
                    "WHERE pseudo_prepid = {} and pipeline_step_id = 11"
                    ).format(pseudo_prepid)
    curs.execute(timesRunQuery)
    timesRun = curs.fetchone()

    if timesRun:
        pipelineID_update = ("UPDATE dragen_pipeline_step SET "
                "submit_time = '{}', "
                "finish_time = '{}', "
                "times_ran = times_ran + 1 "
                "WHERE pseudo_prepid = {} and pipeline_step_id = 11"
                ).format(submitTime,finishTime,pseudo_prepid)
        if debug:
           print pipelineID_update
        curs.execute(pipelineID_update)
    else:
        pipelineID_insert = ("INSERT INTO dragen_pipeline_step "
                "(pseudo_prepid,pipeline_step_id,submit_time,finish_time,times_ran,finished) "
                "VALUES ({},11,'{}','{}',1,1)").format(pseudo_prepid,submitTime,finishTime)
        if debug:
            print pipelineID_insert
        curs.execute(pipelineID_insert)

def update_queue(curs,dragen_id,debug):
    insert_query = ("INSERT INTO tmp_dragen "
            "SELECT * FROM dragen_queue WHERE dragen_id={0}"
            ).format(dragen_id)
    rm_query = ("DELETE FROM {0} WHERE dragen_id={1}").format("dragen_queue",dragen_id)

    curs.execute(insert_query)
    curs.execute(rm_query)
    if debug:
        print insert_query
        print rm_query

def set_seqtime(curs,sample):
    query = ("SELECT FROM_UNIXTIME(Seqtime) "
            "FROM Flowcell WHERE FCIllumID='{first_flowcell}'").format(**sample.metadata)
    curs.execute(query)
    seqtime = curs.fetchone()

    if seqtime: #In case there is not flowcell information
        seqtime = seqtime[0].date().isoformat() #ISO8601 format for RGDT field
    else:
        seqtime = '1970-1-1'
    sample.set('seqtime',seqtime)

def single_sample_setup(curs,sample,parameters,debug):

    setup_dir(curs,sample,debug)
    create_align_config(sample)

def setup_dir(curs,sample,debug):
    mkdir_cmd = ['mkdir','-p',sample.metadata['script_dir']]
    subprocess.call(mkdir_cmd)
    mkdir_cmd = ['mkdir','-p',sample.metadata['log_dir']]
    subprocess.call(mkdir_cmd)
    mkdir_cmd = ['mkdir','-p',sample.metadata['fastq_dir']]
    subprocess.call(mkdir_cmd)

    """Removes any fastq.gz files in the fastq folder. Symlink only fastqs
        should be in this folder.  This is just in case the sample was run
        through once, created symlinks but later was sequenced again.  All
        output is hidden (written to /dev/null/"""
    dev_null = open(os.devnull, 'w')
    files = glob(str(sample.metadata['fastq_dir']) + '/*fastq.gz')
    for file in files:
        #remove any symlinked fastq.gz files
        os.remove(file)

    first_read1 = get_first_read(sample,1,debug)
    first_read2 = get_first_read(sample,2,debug)
    fastq_counter=0
    for fastq_loc in sample.metadata['fastq_loc']:
        fastqs = glob(fastq_loc + '/*_R1_*fastq.gz')
        for fastq in fastqs:
            if 'fastq16-rsync' in fastq:
                pass
            else:
                fastq_counter+=1
                new_fastq_read1 = first_read1.split('/')[-1].replace(
                        '001.fastq.gz','{0:03d}.fastq.gz'.format(fastq_counter))
                new_fastq_read2 = first_read2.split('/')[-1].replace(
                        '001.fastq.gz','{0:03d}.fastq.gz'.format(fastq_counter))
                ln_cmd1 = ['ln','-s',fastq,sample.metadata['fastq_dir']+'/'+new_fastq_read1]
                ln_cmd2 = ['ln','-s',fastq.replace('_R1_','_R2_'),sample.metadata['fastq_dir']+'/'+new_fastq_read2]

                if fastq_counter == 1:
                    sample.set('first_fastq1',sample.metadata['fastq_dir']+'/'+new_fastq_read1)
                    sample.set('first_fastq2',sample.metadata['fastq_dir']+'/'+new_fastq_read2)
                    sample.set('first_lane',sample.metadata['lane'][0][0][0])
                    sample.set('first_flowcell',sample.metadata['lane'][0][0][1])

                subprocess.call(ln_cmd1)
                subprocess.call(ln_cmd2)
    check_Fastq_Total_Size(sample,debug)
    set_seqtime(curs,sample)

def check_Fastq_Total_Size(sample,debug):
    fastq_dir = sample.metadata['fastq_dir']
    sample_type = sample.metadata['sample_type']
    fastq_filesize_sum = 0
    for fastq in glob(fastq_dir + '/*gz'):
        fastq_filesize = os.stat(os.path.realpath(fastq)).st_size
        fastq_filesize_sum += fastq_filesize

    if sample_type == 'genome':
        if fastq_filesize_sum < 53687091200: # < 50GB
            print fastq_filesize_sum
            raise Exception, "Sum of fastq files is too small for a {} sample!".format(sample_type)
    elif sample_type == 'exome':
        if fastq_filesize_sum > 32212254720: # > 30GB
            print fastq_filesize_sum
            raise Exception, "Sum of fastq files is too big for a {} sample!".format(sample_type)
    elif sample_type == 'rnaseq':
        if fastq_filesize_sum > 32212254720: # > 30GB
            print fastq_filesize_sum
            raise Exception, "Sum of fastq files is too big for a {} sample!".format(sample_type)
    elif sample_type == 'custom_capture':
        if fastq_filesize_sum > 10737418240: # > 10GB
            print fastq_filesize_sum
            raise Exception, "Sum of fastq files is too big for a {} sample!".format(sample_type)
    else:
        raise Exception, 'Unhandled sample_type found: {}!'.format(sample_type)
    if debug:
        print 'Fastq size',fastq_filesize_sum,float(fastq_filesize_sum)/1024/1024/1024

def get_first_read(sample,read_number,debug):
    #Using first fastq as template for all fastq.gz
    if debug:
        print sample.metadata['fastq_loc']
    first_fastq_loc = sample.metadata['fastq_loc'][0]
    if debug:
        print '{0}/*L00*_R{1}_001.fastq.gz'.format(first_fastq_loc,read_number)
    read = glob('{0}/*L00*_R{1}_001.fastq.gz'.format(first_fastq_loc,read_number))
    #The fastqs might still be on the quantum tapes
    if read == []:
        raise Exception, "Fastq file not found!"
    else:
        return sorted(read)[0]

def get_next_sample(curs,debug):
    #remove any symlinked fastq.gz files 
    query = ("SELECT sample_name,sample_type,pseudo_prepid,capture_kit,dragen_id "
        "FROM dragen_queue "
        "WHERE PRIORITY < 99 "
        "ORDER BY PRIORITY ASC LIMIT 1 ")

    curs.execute(query)
    dragen_id = curs.fetchone()
    if debug:
        print query
        print 'Dragen_queue info: {0}'.format(dragen_id)
    return dragen_id

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-s", "--sample", action="append", nargs="+",
                        help="Specify sample")
    group.add_argument("--samples", type=open,
                        help="Provide a file with a list of samples with one "
                        "sample per line")
    group.add_argument("-a", "--auto", default=False, action="store_true",
                        help="Run pipeline in automated mode")
    parser.add_argument("-d", "--debug", default=False, action="store_true",
                        help="Verbose output for debugging")
    # --> Add code to accept regexp in auto mode <--#
    # --> Add code to for custom config files <--#
    parser.add_argument("--dontexecute", default=False, action="store_true",
                        help="Perform setup but do not start a Dragen run")
    parser.add_argument("--test", default=False, action="store_true",
                        help="Query and updates to the database occur on the "
                        "test server")
    parser.add_argument('--version', '-v', action='version',
                        version='%(prog)s v1.0')
    args=parser.parse_args()

    if args.sample:
        samples = args.sample
    elif args.samples:
        samples = args.samples
    else:
        samples = args.auto

    if args.test:
        args.database="testDB"
    else:
        args.database="sequenceDB"


    main(samples, args.debug, args.dontexecute, args.database)

