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
#from create_config import create_config
from create_align_config import create_align_config
#from create_gvcf_config import create_gvcf_config
#from create_joint_call_config import create_joint_call_config
from ConfigParser import SafeConfigParser
from datetime import datetime
from dragen_sample import dragen_sample
from glob import glob

def main(samples, debug, dontexecute, gvcf, database):
    CNF = "/nfs/goldstein/software/dragen/dragen.cnf" # defaults file for pipeline
    config_parser = SafeConfigParser()
    config_parser.read(CNF)
    parameters = {}
    for parameter in ("PYTHON", "REFERENCE_GENOME", "JAVA_1_8", "PIGZ",
            "DB_CONFIG_FILE", "PICARD_JAR", "DRAGEN_CONFIG"):
        parameters[parameter] = config_parser.get("dragen_pipe", parameter.upper())

    if debug:
        print "Parameters used: {0}".format(parameters)
        if gvcf:
            print 'gVCF mode activated'

    db = MySQLdb.connect(db=database, read_default_group="clientsequencedb",
            read_default_file=parameters["DB_CONFIG_FILE"])
    curs = db.cursor()

    #Automated run
    if samples == True:
        info = get_next_sample(curs,debug)
        if info is None:
            print "No samples were found"
            sys.exit()

        while info is not None:
            dragen_id = info[4]
            sample = dragen_sample(info[0],info[1],info[2],info[3],curs)
            single_sample_setup(curs,sample,parameters,gvcf,debug)
            run_sample(curs,sample,dragen_id,dontexecute,debug)
            info = get_next_sample(curs,debug)

    #Run through a sample file
    else:
        for line in samples.readlines():
            info = line.strip().split('\t')
            prep_query = ("SELECT pseudo_prepid "
                "FROM seqdbClone "
                "WHERE CHGVID='{0}' AND seqtype ='{1}' "
                "ORDER BY prepid desc limit 1"
                ).format(info[0],info[1])
            curs.execute(prep_query)
            pseudo_prepid = curs.fetchone()
            sample = dragen_sample(info[0],info[1],pseudo_prepid[0],info[2],curs)
            single_sample_setup(curs,sample,parameters,gvcf,debug)
            run_sample(curs,sample,0,dontexecute,debug)

def run_sample(curs,sample,dragen_id,dontexecute,debug):
    cmd = ['dragen', '-f', '-v', '-c', sample.metadata['conf_file']]
    if debug:
        print ' '.join(cmd)

    dragen_stderr = open(sample.metadata['dragen_stderr'],'a')
    if dontexecute == False:
        if dragen_id != 0: # All manually run samples will have a dragen_id of 0
            update_queue(curs,dragen_id,debug)

        with open(sample.metadata['dragen_stdout'],'a') as dragen_stdout:

            process = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=dragen_stderr)
            for line in iter(process.stdout.readline, ''):
                if debug:
                    sys.stdout.write(line)
                dragen_stdout.write(line)

            process.communicate()
            error_code = process.wait()

        if debug:
           print "Dragen error code: {0}".format(error_code)
        if error_code == 0:
            rm_query = "DELETE FROM {0} WHERE dragen_id={1}".format("tmp_dragen",dragen_id)
            if debug:
                print rm_query
            curs.execute(rm_query)
        else:
            print "Sample was not deleted from tmp_dragen due to error code: {0}".format(error_code)

        dragen_stdout.close()
        dragen_stderr.close()
    else:
        print "Sample {sample_name} config complete.".format(**sample.metadata)

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

def single_sample_setup(curs,sample,parameters,gvcf,debug):

    setup_dir(curs,sample,debug)

    create_align_config(sample)
    #create_config(sample,gvcf)
    #create_gvcf_config(sample)
    #create_joint_call_config(sample)
    #create_post_dragen_shell(sample,parsed_CNF)

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
        os.remove(file)

    first_read1 = get_first_read(sample,1,debug)
    first_read2 = get_first_read(sample,2,debug)
    fastq_counter=0
    for fastq_loc in sample.metadata['fastq_loc']:
        fastqs = glob(fastq_loc + '/*R1*fastq.gz')
        for fastq in fastqs:
            fastq_counter+=1
            new_fastq_read1 = first_read1.split('/')[-1].replace(
                    '001.fastq.gz','{0:03d}.fastq.gz'.format(fastq_counter))
            new_fastq_read2 = first_read2.split('/')[-1].replace(
                    '001.fastq.gz','{0:03d}.fastq.gz'.format(fastq_counter))
            ln_cmd1 = ['ln','-s',fastq,sample.metadata['fastq_dir']+'/'+new_fastq_read1]
            ln_cmd2 = ['ln','-s',fastq.replace('R1','R2'),sample.metadata['fastq_dir']+'/'+new_fastq_read2]

            if fastq_counter == 1:
                sample.set('first_fastq1',sample.metadata['fastq_dir']+'/'+new_fastq_read1)
                sample.set('first_fastq2',sample.metadata['fastq_dir']+'/'+new_fastq_read2)
                sample.set('first_lane',sample.metadata['lane'][0][0][0])
                sample.set('first_flowcell',sample.metadata['lane'][0][0][1])

            subprocess.call(ln_cmd1)
            subprocess.call(ln_cmd2)
    set_seqtime(curs,sample)

def get_first_read(sample,read_number,debug):
    #Using first fastq as template for all fastq.gz
    first_fastq_loc = sample.metadata['fastq_loc'][0]
    if debug:
        print '{0}/*L00{1}_R{2}_001.fastq.gz'.format(first_fastq_loc,sample.metadata['lane'][0][0][0],read_number)
    read = glob('{0}/*L00{1}_R{2}_001.fastq.gz'.format(first_fastq_loc,sample.metadata['lane'][0][0][0],read_number))
    #The fastqs might still be on the quantum tapes
    if read == []:
        read = glob(('/stornext/seqfinal/casava1.8/whole_exome/{0}/{1}/*L00{2}_R{3}_001.fastq.gz'
            ).format(sample.metadata['sample_name'],sample.metadata['lane'][0][0][1],
                sample.metadata['lane'][0][0][0],read_number))

        if read == []:
            print sample.metadata['sample_name']

            raise Exception, "Fastq file not found!"
        else:
            """fastq_loc was based off of database so needs to be set to the
                quantum location"""
            sample.set('fastq_loc',glob(('/stornext/seqfinal/casava1.8/whole_exome/{0}/*XX'
                ).format(sample.metadata['sample_name'])))
    return read[0]

def get_next_sample(curs,debug):
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
    parser.add_argument("--gvcf", default=False, action="store_true",
                        help="Run dragen and product a gVCF")
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


    main(samples, args.debug, args.dontexecute, args.gvcf, args.database)

