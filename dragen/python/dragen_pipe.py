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
from create_config import create_config
#from create_gvcf_config import create_gvcf_config
#from create_joint_call_config import create_joint_call_config
from ConfigParser import SafeConfigParser
from datetime import datetime
from dragen_sample import dragen_sample
from glob import glob

def main(samples, debug, execute, database):
    CNF = "/nfs/goldstein/software/dragen/dragen.cnf" # defaults file for pipeline
    config_parser = SafeConfigParser()
    config_parser.read(CNF)
    parameters = {}
    for parameter in ("PYTHON", "REFERENCE_GENOME", "JAVA_1_8", "PIGZ",
            "DB_CONFIG_FILE", "PICARD_JAR", "DRAGEN_CONFIG"):
        parameters[parameter] = config_parser.get("dragen_pipe", parameter.upper())

    if debug:
        print "Parameters used: {0}".format(parameters)

    db = MySQLdb.connect(db=database, read_default_group="clientsequencedb",
            read_default_file=parameters["DB_CONFIG_FILE"])
    curs = db.cursor()

    if samples == True: #Automated run

        info = get_next_sample(curs,debug)
        if info is None:
            print "No samples were found"
            sys.exit()

        while info is not None:
            dragen_id = info[4]
            sample = dragen_sample(info[0],info[1],info[2],info[3],curs)
            single_sample_setup(curs,sample,parameters)
            run_sample(sample)
            query = "DELETE FROM dragen_queue WHERE dragen_id={0}".format(dragen_id)
            if debug:
                print query
            curs.execute(query)
            info = get_next_sample(curs,debug)

        pass
    else:
        pass
        #dragen_id = get_sample(sample_name,sample_type,capture_kit)
        #single_sample_setup(dragen_id)

def run_sample(sample):
    cmd = ['dragen', '-f', '-v', '-c', sample.metadata['conf_file'], sample.metadata['out_file']]
    subprocess.call(cmd, stdout=sample.metadata['dragen_stdout'],stderr=sample.metadata['dragen_stderr'])

def single_sample_setup(curs,sample,parameters):

    setup_dir(sample)
    query = ("SELECT FROM_UNIXTIME(Seqtime) "
            "FROM Flowcell WHERE FCIllumID='{first_flowcell}'"
            ).format(**sample.metadata)
    curs.execute(query)
    seqtime = curs.fetchone()
    sample.set('seqtime',seqtime)

    create_config(sample)
    #create_gvcf_config(sample)
    #create_joint_call_config(sample)
    #create_post_dragen_shell(sample,parsed_CNF)

def setup_dir(sample):
    mkdir_cmd = ['mkdir','-p',sample.metadata['script_dir']]
    subprocess.call(mkdir_cmd)
    mkdir_cmd = ['mkdir','-p',sample.metadata['log_dir']]
    subprocess.call(mkdir_cmd)
    mkdir_cmd = ['mkdir','-p',sample.metadata['fastq_dir']]
    subprocess.call(mkdir_cmd)

    first_read1 = get_first_read(sample,1)
    first_read2 = get_first_read(sample,2)
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
            ln_cmd2 = ['ln','-s',fastq,sample.metadata['fastq_dir']+'/'+new_fastq_read2]

            if fastq_counter == 1:
                sample.set('first_fastq1',sample.metadata['fastq_dir']+'/'+new_fastq_read1)
                sample.set('first_fastq2',sample.metadata['fastq_dir']+'/'+new_fastq_read2)
                sample.set('first_lane',sample.metadata['lane'][0][0][0])
                sample.set('first_flowcell',sample.metadata['lane'][0][0][1])

            #print ln_cmd1
            #print ln_cmd2
            subprocess.call(ln_cmd1)
            subprocess.call(ln_cmd2)


def get_first_read(sample,read_number):
    #Using first fastq as template for all fastq.gz 
    first_fastq_loc = sample.metadata['fastq_loc'][0]
    read = glob('{0}/*L00{1}_R{2}_001.fastq.gz'.format(first_fastq_loc,sample.metadata['lane'][0][0][0],read_number))
    return read[0]


def get_next_sample(curs,debug):
    query = ("SELECT sample_name,sample_type,pseudo_prepid,capture_kit,dragen_id "
        "FROM dragen_queue "
        "ORDER BY PRIORITY DESC LIMIT 1 ")
    if debug:
        print query

    curs.execute(query)
    dragen_id = curs.fetchone()
    return dragen_id

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-s", "--sample", action="append", nargs="+",
                        help="Specify sample")
    group.add_argument("--samples", type=open,
                        help="Provide a file with a list of samples with one "
                        "sample per line")
    group.add_argument("-p", "--prepid", action="append", nargs="+",
                        help="Specify prepid")
    group.add_argument("--prepids", type=open,
                        help="Provide a file with a list of prepid with one "
                        "prepid per line")
    group.add_argument("-a", "--auto", default=False, action="store_true",
                        help="Run pipeline in automated mode")
    parser.add_argument("-d", "--debug", default=False, action="store_true",
                        help="Verbose output for debugging")
    # --> Add code to accept regexp in auto mode <--#
    # --> Add code to for custom config files <--#
    parser.add_argument("--execute", default=True, action="store_false",
                        help="Automatically run dragen")
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
    elif args.prepid:
        samples = args.prepid
    elif args.prepids:
        samples = args.prepids
    else:
        samples = args.auto

    if args.test:
        args.database="testDB"
    else:
        args.database="sequenceDB"


    main(samples, args.debug, args.execute, args.database)

