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
from ConfigParser import SafeConfigParser
from create_align_config import create_align_config
from datetime import datetime
from dragen_sample import dragen_sample
from glob import glob

CNF = "/nfs/goldstein/software/dragen/dragen.cnf" # defaults file for pipeline
config_parser = SafeConfigParser()
config_parser.read(CNF)
parameters = {}
for parameter in ("DB_CONFIG_FILE","DRAGEN_CONFIG"):
    parameters[parameter] = config_parser.get("dragen_pipe", parameter.upper())

def main(samples, debug, dontexecute, seqscratch_drive):
    if debug:
        print "Parameters used: {0}".format(parameters)

    #Automated run
    sample_name, sample_type, pseudo_prepid, capture_kit, dragen_id = get_next_sample(debug)

    while sample_name is not None:
        sample = dragen_sample(sample_name,sample_type,pseudo_prepid,capture_kit,get_curs()[0])
        sample.metadata['output_dir'] = ('/nfs/{}/ALIGNMENT/BUILD37/DRAGEN/{}/{}.{}/'
                                        ).format(seqscratch_drive,sample_type.upper(),
                                                 sample_name,pseudo_prepid)
        single_sample_setup(sample,parameters,debug)
        run_sample(sample,dragen_id,dontexecute,parameters,debug)

        try:
            sample_name, sample_type, pseudo_prepid, capture_kit, dragen_id = get_next_sample(debug)
        except:
            print "No more samples in the queue"
            sys.exit()

def get_curs():
    db = MySQLdb.connect(db=database, read_default_group="clientsequencedb",
        read_default_file=parameters["DB_CONFIG_FILE"])
    curs = db.cursor()
    return curs,db

def run_query(query):
    curs,db = get_curs()
    curs.execute(query)
    results = curs.fetchall()
    db.commit()
    db.close()
    return results

def run_sample(sample,dragen_id,dontexecute,parameters,debug):
    output_dir = sample.metadata['output_dir']
    pseudo_prepid = sample.metadata['pseudo_prepid']
    existing_bams = glob("{}/*.sam".format(output_dir))
    submitTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    if existing_bams == []:
        update_pipeline_step_id(sample,submitTime,'Null','started',debug)
        dragenCmd = ['dragen','-f','-v','-c',sample.metadata['conf_file'],
               '--watchdog-active-timeout','7000']
        if debug:
            print ' '.join(dragenCmd)

        if dragen_id != 0: # All manually run samples will have a dragen_id of 0
            update_queue(dragen_id,debug)

        dragen_stderr = open(sample.metadata['dragen_stderr'],'a')
        with open(sample.metadata['dragen_stdout'],'a') as dragen_stdout:
            process = subprocess.Popen(dragenCmd,stdout=subprocess.PIPE,stderr=dragen_stderr)
            for line in iter(process.stdout.readline, ''):
                if debug:
                    sys.stdout.write(line)
                dragen_stdout.write(line)
            process.communicate()
            error_code = process.wait()
            subprocess.call(['chmod','-R','775','{}'.format(output_dir)])
        if error_code == 0:
            finishTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_query = ("INSERT INTO gatk_queue "
                            "SELECT * FROM tmp_dragen WHERE dragen_id={0}"
                           ).format(dragen_id)
            rm_query = "DELETE FROM {0} WHERE dragen_id={1}".format("tmp_dragen",dragen_id)

            if debug:
                print insert_query
                print rm_query
            run_query(insert_query)
            run_query(rm_query)
            update_pipeline_step_id(sample,submitTime,finishTime,'completed',debug)
            update_dragen_metadata(sample,debug)

        else:
            print ("Sample was not deleted from tmp_dragen due to "
                   "error code: {0}").format(error_code)

        dragen_stdout.close()
        dragen_stderr.close()
    else:
        print ("Sample {sample_name} bam file already exists!"
              ).format(sample_name=sample.metadata['sample_name'])
        sys.exit()
    update_queue(dragen_id,debug)

def update_dragen_metadata(sample,debug):
    seqscratch_drive = sample.metadata['output_dir'].split('/')[2]

    query = ("SELECT * from dragen_sample_metadata "
            "WHERE pseudo_prepid = {pseudo_prepid} "
            ).format(**sample.metadata)
    dbInfo = run_query(query)
    if dbInfo:
        pass
    else:
        insertQuery = ("INSERT INTO dragen_sample_metadata "
                       "(sample_name,pseudo_prepid,sample_type,"
                       "capture_kit,priority,seqscratch_drive) "
                       "VALUES ('{sample_name}',{pseudo_prepid},"
                       "'{sample_type}','{capture_kit}',{priority},'{seqscratch_drive}') "
                      ).format(seqscratch_drive=seqscratch_drive,**sample.metadata)
        print insertQuery
        run_query(insertQuery)

def get_pipeline_version():
    version = subprocess.check_output(["/nfs/goldstein/software/git-2.5.0/bin/git", "describe", "--tags"]).strip()
    if version:
        return version
    else:
            raise ValueError("Could not get the version # of the pipeline; "
                             "maybe run it from a directory in the repo?")

def perform_DPS_update(version,submitTime,finishTime,timesRan,
                       status,pseudo_prepid,alignment_step_id,debug):
    pipelineID_update = ("UPDATE dragen_pipeline_step SET "
                         "version = '{}', "
                         "submit_time = '{}', "
                         "finish_time = '{}', "
                         "times_ran ={}, "
                         "step_status = '{}' "
                         "WHERE pseudo_prepid = {} "
                         "and pipeline_step_id = {}"
                        ).format(version,submitTime,finishTime,
                                 timesRan,status,pseudo_prepid,
                                 alignment_step_id)
    if debug:
       print pipelineID_update
    run_query(pipelineID_update)

def perform_DPS_insert(version,submitTime,finishTime,timesRan,
                       status,pseudo_prepid,alignment_step_id,debug):
        pipelineID_insert = ("INSERT INTO dragen_pipeline_step "
                             "(version,pseudo_prepid,pipeline_step_id,"
                             "submit_time,finish_time,times_ran,step_status) "
                             "VALUES ('{}',{},{},'{}','{}',1,'{}')"
                            ).format(version,pseudo_prepid,
                                     alignment_step_id,submitTime,
                                     finishTime,status)
        if debug:
            print pipelineID_insert
        run_query(pipelineID_insert)

def update_pipeline_step_id(sample,submitTime,finishTime,status,debug):
    version = get_pipeline_version()
    alignment_step_query = ("SELECT id from dragen_pipeline_step_desc "
                            "WHERE step_name='DragenAlignment'")
    alignment_step_id = run_query(alignment_step_query)[0][0]
    pseudo_prepid = sample.metadata['pseudo_prepid']
    timesRanQuery = ("SELECT times_ran FROM dragen_pipeline_step "
                     "WHERE pseudo_prepid = {} and pipeline_step_id = {}"
                    ).format(pseudo_prepid,alignment_step_id)
    timesRan = run_query(timesRanQuery)
    print timesRanQuery,timesRan
    if timesRan and status == 'started':
        timesRan = int(timesRan[0][0]) +1
        perform_DPS_update(version,submitTime,finishTime,timesRan,
                           status,pseudo_prepid,alignment_step_id,debug)

    elif timesRan and status == 'completed':
        timesRan = int(timesRan[0][0])
        perform_DPS_update(version,submitTime,finishTime,timesRan,
                           status,pseudo_prepid,alignment_step_id,debug)

    else:
        perform_DPS_insert(version,submitTime,finishTime,timesRan,
                           status,pseudo_prepid,alignment_step_id,debug)

def update_queue(dragen_id,debug):
    insert_query = ("INSERT INTO tmp_dragen "
                    "SELECT * FROM dragen_queue WHERE dragen_id={0}"
                   ).format(dragen_id)
    rm_query = ("DELETE FROM {0} WHERE dragen_id={1}").format("dragen_queue",dragen_id)
    run_query(insert_query)
    run_query(rm_query)
    if debug:
        print insert_query
        print rm_query

def set_seqtime(sample):
    query = ("SELECT FROM_UNIXTIME(Seqtime) "
             "FROM Flowcell WHERE FCIllumID='{first_flowcell}'"
            ).format(**sample.metadata)
    seqtime = run_query(query)

    if seqtime:
        seqtime = seqtime[0][0].date().isoformat() #ISO8601 format for RGDT field
    else: #In case there is not flowcell information (Ex. old and external samples)
        seqtime = '1970-1-1'
    sample.set('seqtime',seqtime)

def single_sample_setup(sample,parameters,debug):
    setup_dir(sample,debug)
    create_align_config(sample)

def setup_dir(sample,debug):
    mkdir_cmd = ['mkdir','-p',sample.metadata['script_dir']]
    subprocess.call(mkdir_cmd)
    mkdir_cmd = ['mkdir','-p',sample.metadata['log_dir']]
    subprocess.call(mkdir_cmd)
    mkdir_cmd = ['mkdir','-p',sample.metadata['fastq_dir']]
    subprocess.call(mkdir_cmd)

    """Removes any fastq.gz files in the fastq folder. Symlink-only fastqs
        should be in this folder.  This is just in case the sample was run
        through once"""
    files = glob(str(sample.metadata['fastq_dir']) + '/*fastq.gz')
    for file in files:
        #remove any symlinked fastq.gz files
        os.remove(file)

    first_read1 = get_reads(sample,1,debug)
    first_read2 = get_reads(sample,2,debug)
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
    set_seqtime(sample)

def check_Fastq_Total_Size(sample,debug):
    fastq_dir = sample.metadata['fastq_dir']
    sample_type = sample.metadata['sample_type']
    fastq_filesize_sum = 0
    for fastq in glob(fastq_dir + '/*gz'):
        fastq_filesize = os.stat(os.path.realpath(fastq)).st_size
        fastq_filesize_sum += fastq_filesize
    print "Sum of Fastq size: {}".format(fastq_filesize_sum)
    if sample_type == 'genome':
        if fastq_filesize_sum < 53687091200: # < 50GB
            userInput = raw_input('Sum of fastq files sizes are too small.  Is this ok? (y)es or (n)o ').lower()
            if userInput.lower() == 'n':
                raise Exception, "Sum of fastq files sizes is too small for a {} sample!".format(sample_type)
    elif sample_type == 'exome':
        if fastq_filesize_sum > 32212254720: # > 30GB
            userInput = raw_input('Sum of fastq files sizes are too big.  Is this ok? (y)es or (n)o ').lower()
            if userInput.lower() == 'n':
                 raise Exception, "Sum of fastq files is too big for a {} sample!".format(sample_type)
        elif fastq_filesize_sum < 1073741824: # < 1GB
            userInput = raw_input('Sum of fastq files sizes are too small.  Is this ok? (y)es or (n)o ').lower()
            if userInput.lower() == 'n':
                 raise Exception, "Sum of fastq files is too small for a {} sample!".format(sample_type)
    elif sample_type == 'rnaseq':
        if fastq_filesize_sum > 32212254720: # > 30GB
            userInput = raw_input('Sum of fastq files sizes are too big.  Is this ok? (y)es or (n)o ').lower()
            if userInput.lower() == 'n':
                 raise Exception, "Sum of fastq files sizes is too big for a {} sample!".format(sample_type)
        elif fastq_filesize_sum < 1073741824: # < 1GB
            userInput = raw_input('Sum of fastq files sizes are too small.  Is this ok? (y)es or (n)o ').lower()
            if userInput.lower() == 'n':
                 raise Exception, "Sum of fastq files is too small for a {} sample!".format(sample_type)
    elif sample_type == 'custom_capture':
        if fastq_filesize_sum > 10737418240: # > 10GB
            userInput = raw_input('Sum of fastq files sizes are too big.  Is this ok? (y)es or (n)o ').lower()
            if userInput.lower() == 'n':
                 raise Exception, "Sum of fastq files sizes is too big for a {} sample!".format(sample_type)
    else:
        raise Exception, 'Unhandled sample_type found: {}!'.format(sample_type)
    if debug:
        print 'Fastq size',fastq_filesize_sum,float(fastq_filesize_sum)/1024/1024/1024

def get_reads(sample,read_number,debug):
    if debug:
        print sample.metadata['fastq_loc']
    fastq_loc = sample.metadata['fastq_loc'][0]
    if debug:
        print '{0}/*L00*_R{1}_001.fastq.gz'.format(fastq_loc,read_number)
    read = glob('{0}/*L00*_R{1}_001.fastq.gz'.format(fastq_loc,read_number))
    #The fastqs might still be on the quantum tapes
    if read == []:
        raise Exception, "Fastq file not found!"
    else:
        return sorted(read)[0]

def get_next_sample(debug):
    #remove any symlinked fastq.gz files 
    query = ("SELECT sample_name,sample_type,pseudo_prepid,capture_kit,dragen_id "
        "FROM dragen_queue "
        "WHERE PRIORITY < 99 "
        "ORDER BY PRIORITY ASC LIMIT 1 ")
    sample_info = run_query(query)
    if debug:
        print query
        print 'Dragen_queue info: {0}'.format(sample_info)
    return sample_info[0]

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
    parser.add_argument("--seqscratch_drive", default='seqscratch_ssd', action="store_true",
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
    global database
    database = args.database
    main(samples, args.debug, args.dontexecute, args.seqscratch_drive)

