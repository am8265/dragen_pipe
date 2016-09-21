#!/nfs/goldstein/software/python2.7.7/bin/python2.7

import argparse
import MySQLdb
import os
import shlex
import string
import sys
import subprocess
import time
from ConfigParser import SafeConfigParser
from random import randint


def main(database,debug):
    CNF = "/nfs/goldstein/software/dragen/dragen.cnf" # defaults file for pipeline
    config_parser = SafeConfigParser()
    config_parser.read(CNF)
    parameters = {}
    for parameter in ("PYTHON", "REFERENCE_GENOME", "JAVA_1_8", "PIGZ",
            "DB_CONFIG_FILE", "PICARD_JAR", "DRAGEN_CONFIG"):
        parameters[parameter] = config_parser.get("dragen_pipe", parameter.upper())
    db = MySQLdb.connect(db=database, read_default_group="clientsequencedb",
            read_default_file=parameters["DB_CONFIG_FILE"])
    curs = db.cursor()

    run_samples(curs,debug)

def run_samples(curs,debug):
    dragen_sample = get_next_sample(curs,debug)
    number_of_jobs = get_avail_cores()
    while dragen_sample and number_of_jobs > 0:
        number_of_jobs = get_avail_cores()

        if debug:
            print "Number of Jobs: {0}".format(number_of_jobs)
            print dragen_sample
        if dragen_sample is None:
            print "No samples were found"
            sys.exit()
        add_to_tmp_gatk(curs,dragen_sample,debug)
        run_sample(curs,dragen_sample,debug)

        number_of_jobs = 0

def add_to_tmp_gatk(curs,dragen_sample,debug):
    dragen_id = dragen_sample[4]
    insert_query = ("INSERT INTO tmp_gatk "
                    "SELECT * FROM gatk_queue WHERE dragen_id={0}"
                    ).format(dragen_id)
    rm_query = "DELETE FROM gatk_queue WHERE dragen_id={}".format(dragen_id)
    if debug:
        print insert_query
        print rm_query

    curs.execute(insert_query)
    curs.execute(rm_query)

def run_sample(curs,dragen_sample,debug):
    python2_7_7_loc = "/nfs/goldstein/software/python2.7.7/bin"
    end_module = "ArchiveSample"
    sample_name = dragen_sample[0]
    sample_type = dragen_sample[1].lower()
    pseudo_prepid = dragen_sample[2]
    capture_kit_info = get_capture_kit_info(curs,dragen_sample,debug)
    
    
    if sample_type == 'custom capture':
        sample_type = sample_type.replace(" ","_")
    export_path_cmd = "export PATH={python2_7_7_loc}:$PATH".format(python2_7_7_loc=python2_7_7_loc)
    export_luigi_config_cmd = 'export LUIGI_CONFIG_PATH=/home/jb3816/github/dragen_pipe/dragen/python/luigi.cfg'
    luigi_cmd = ( "{python2_7_7_loc}/python2.7 "
    "-m luigi --module gatk_pipe {end_module} "
    "--sample-name {sample_name} "
    "--capture-kit-bed {capture_kit_bed} "
    "--sample-type {sample_type} "
    "--pseudo-prepid {pseudo_prepid} "
    "--no-tarball "
    "--poll-time 60 "
    "--worker-wait-interval 60 "
    "--workers 3 --scheduler-remove-delay 86400"
    ).format(capture_kit_bed=capture_kit_info[0],
            end_module=end_module,
            sample_type=sample_type,
            sample_name=sample_name,
            python2_7_7_loc=python2_7_7_loc,
            pseudo_prepid=pseudo_prepid
            )
    if debug:
        print luigi_cmd
    print export_path_cmd
    #subprocess.check_call(shlex.split(export_path_cmd),shell=True)
    #subprocess.check_call(shlex.split(export_luigi_config_cmd),shell=True)
    proc = subprocess.Popen(shlex.split(luigi_cmd),shell=True)
    return proc

def get_capture_kit_info(curs,dragen_sample,debug):
    capture_kit = dragen_sample[3]
    if dragen_sample[1] == 'genome':
         query = ("SELECT DISTINCT region_file_lsrc,name "
            "FROM captureKit WHERE prepT_name='Roche SeqCap EZ V3'")
    else:
        query = ("SELECT DISTINCT region_file_lsrc,name "
            "FROM captureKit WHERE prepT_name='{capture_kit}'"
            ).format(capture_kit=capture_kit)

    if debug:
        print query
    curs.execute(query)
    capture_kit_bed = curs.fetchone()

    #changed to list since the script will alter this value if its not custom capture
    return list(capture_kit_bed)

def get_next_sample(curs,debug):
    query = ("SELECT sample_name,sample_type,pseudo_prepid,capture_kit,dragen_id "
        "FROM gatk_queue "
        "WHERE PRIORITY < 99 "
        "ORDER BY PRIORITY ASC LIMIT 1 ")

    curs.execute(query)
    dragen_sample = curs.fetchone()
    if debug:
        print query
        print 'Dragen_queue ID: {0}'.format(dragen_sample[4])
    return dragen_sample


def get_avail_cores():
    cores_left_avail = 200
    cores_per_pipeline = 4

    qstat_cmd = ['qstat','-g','c']
    output = subprocess.check_output(qstat_cmd).splitlines()[2]
    total_avail_cores = int([s for s in output.split(' ') if s][5])

    if total_avail_cores > cores_left_avail:
        avail_cores = total_avail_cores - cores_left_avail
        return avail_cores/cores_per_pipeline
    else:
        raise Exception, "Not enough cores are available"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-d", "--debug", default=False, action="store_true",
                        help="Verbose output for debugging")
    parser.add_argument("--dontexecute", default=False, action="store_true",
                        help="Perform setup but do not start a Dragen run")
    parser.add_argument("--test", default=False, action="store_true",
                        help="Query and updates to the database occur on the "
                        "test server")
    parser.add_argument('--version', '-v', action='version',
                        version='%(prog)s v1.0')
    args=parser.parse_args()
    if args.test:
        args.database="testDB"
    else:
        args.database="sequenceDB"

    main(args.database,args.debug)
