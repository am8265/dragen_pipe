#!/nfs/goldstein/software/python2.7.7/bin/python2.7

import argparse
import MySQLdb
import os
import shlex
import string
import sys
import subprocess
import time
import random
from dragen_globals import get_connection
from datetime import datetime

##############################################################################
####     This script will submit a luigi job to run the gatk pipeline     ####
####     qc metrics and update to database steps for a single sample      ####
####              Add this to your crontab for automation                 ####
##############################################################################                                                                                                                        
def execute_query(database,query,fetch):
    """ Execute db fetch statements

    database : The database to use for the query
    query : The query to be performed
    fetch : Bool; Indicates whether to fetch the results incase of a SELECT statement
    
    return : Nothing
    """

    tries = 3
    for i in range(tries):
        try:
            db = get_connection(database)
            cur = db.cursor()
            cur.execute(query)
            if fetch == True: ## For select statements 
                result = cur.fetchone()
                return result
            else: ## For insert and delete statements
                db.commit()
        except MySQLdb.Error, e:
            if i == tries - 1:
                raise Exception("Error %d IN CONNECTION: %s" %(e.args[0],e.args[1]))
            else:
                time.sleep(random.randint(20,120)) ## Wait a random duration between 20 secs to a couple of minutes before trying again
                continue 
        finally:
            db.close()
            
def get_sample(args):
    """ Get the next highest priority sample from the gatk_queue table
    """

    query = ("SELECT * "
             "FROM gatk_queue "
             "WHERE PRIORITY < 99 "
             "ORDER BY PRIORITY ASC")
    
    result = execute_query(args.database,query,fetch=True)
    return result

def add_to_tmp(args,sample_info):
    """ Add to the temp queue
    """
    
    dragen_id,pseudo_prepid,sample_name,sample_type,capture_kit,priority,failed = sample_info[0:]
    query = (""" INSERT INTO tmp_gatk ('dragen_id','pseudo_prepid','sample_name','sample_type','capture_kit','priority','failed') VALUES({0},{1}.{2},{3},{4},{5},{6})""".format(dragen_id,pseudo_prepid,sample_name,sample_type,capture_kit,priority,failed))
    execute_query(args.database,query,fetch=False)

def remove_sample(args,sample_info):
    """ Remove the popped sample from the gatk_queue table
    """
    
    dragen_id,pseudo_prepid,sample_name,sample_type,capture_kit,priority,failed = sample_info[0:]
    query = (""" DELETE FROM gatk_queue WHERE sample_name = {0} AND sample_type = {1} AND pseudo_prepid = {2} limit 1""".format(sample_name,sample_type,pseudo_prepid))
    execute_query(args.database,query,fetch=False)

def create_logdir(baselogdir,sample_info):
    """ Create luigi stdout and stderr logging directory for a sample
    
    logdir : string ; Base logging directory
    sample_info : list ; Sample information from gatk queue 

    return : string ; Final path to the sample log directory 
    """

    dragen_id,pseudo_prepid,sample_name,sample_type,capture_kit,priority,failed = sample_info[0:]
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d")
    sample_log_dir = os.path.join(timestamp+'_'+'logs',sample_name+'_'+str(pseudo_prepid))
    final_log_dir = os.path.join(baselogdir,sample_log_dir)
    if not os.path.exists(final_log_dir):
        os.makedirs(final_log_dir)
        print final_log_dir
    return final_log_dir

def get_capturekit_info(sample_type,capture_kit):
    """ Get the capturebed file location
    sample_type : string ; exome/genome/custom capture etc.
    capture_kit : string ; the capture kit type , e.g. Roche SeqCap EZ V3, 50MB, etc.

    returns : string ; the capture bed file location 
    """
    
    if sample_type.upper() == 'GENOME':
         query = ("SELECT region_file_lsrc "
            "FROM captureKit WHERE prepT_name='Roche SeqCap EZ V3' "
            "AND chr = 'all' ")
         result = execute_query(args.database,query,fetch=True)
         return result[0]
    else:
        query = ("SELECT region_file_lsrc "
                 "FROM captureKit WHERE prepT_name='{0}' "
                 "AND chr = 'all'").format(capture_kit)
        result = execute_query(args.database,query,fetch=True)
        return result[0]

def is_complete():
    return False

def run_pipeline(args,sample_info,final_log_dir):
    """ Run a luigi job     
    """
    
    dragen_id,pseudo_prepid,sample_name,sample_type,capture_kit,priority,failed = sample_info[0:]
    pseudo_prepid = str(pseudo_prepid)
    capturekit_bed = get_capturekit_info(sample_type,capture_kit)
    out_file = os.path.join(final_log_dir,sample_name+'_'+pseudo_prepid+'.stdout')
    error_file = os.path.join(final_log_dir,sample_name+'_'+pseudo_prepid+'.stderr')
    cmd = (""" luigi --module gatk_pipe ArchiveSample --sample-name {0} --pseudo-prepid {1} --capture-kit-bed  {2} --sample-type {3} --dont-remove-tmp-dir --no-tarball --poll-time 120 --workers 2 --worker-wait-interval 180 --scheduler-remove-delay 86400 >> {4} 2>>{5}""".format(sample_name,pseudo_prepid,capturekit_bed,sample_type,out_file,error_file))
    print cmd
    
    """
    proc = subprocess.Popen(cmd,shell=True)
    proc.wait()
    if proc.returncode: ## Non zero return code
        ## There are a few options here ;
        ## 1. Re-run the pipeline, luigi already does this and it might be worth increasing the retry attempts
        ## 2. Log the failure in database and rerun at some other point in time, this info can
        ##    already be obtain from the dragen_pipeline_step table
        ## 3. Wait for a random time and rerun the command again, this might be a good option for now
        
        ## Note : Recursion is also a good option, as long as this wrapper is bug free, it will give us
        ## a window of opppurtunity to fix the pipeline code before recursing for another try upon
        ## failure. This would be the final solution for complete automation  
        
        ## Going with option 3 for now, go to sleep :)
        time.sleep(random.randint(3000,10000)) ## Wait a random time between 50 and 166 minutes

        ## We could be on to the next day, so check for the existence and create a logdir if it
        ## does not exist yet for this sample
        final_log_dir = create_logdir(args.baselogdir,sample_info)
        out_file = os.path.join(final_log_dir,sample_name+'_'+pseudo_prepid+'.stdout')
        error_file = os.path.join(final_log_dir,sample_name+'_'+pseudo_prepid+'.stderr')
    """
        #cmd = (""" luigi --module gatk_pipe ArchiveSample --sample-name {0} --pseudo-prepid {1} --capture-kit-bed  {2} --sample-type {3} --dont-remove-tmp-dir --no-tarball --poll-time 120 --workers 2 --worker-wait-interval 180 --scheduler-remove-delay 86400 >> {4} 2>>{5}""".format(sample_name,pseudo_prepid,capturekit_bed,sample_type,out_file,error_file))
        
    
def main(args):
    """ The main function
    """
    
    while True:
        sample_info = get_sample(args)
        #add_to_tmp(args,sample_info)
        final_log_dir = create_logdir(args.baselogdir,sample_info)
        run_pipeline(args,sample_info,final_log_dir)
        #if is_complete(sample_info):
            #remove_sample(args,sample_info)
                           
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Luigi wrapper script',description="This script is a wrapper for running the new dragen variant calling and qc pipeline. It is created for a very specific use case i.e. getting sample information from the table gatk_queue in sequenceDB , adding it to the table tmp_queue , deleting that entry from gatk_queue and submitting a luigi job for the sample. To automate this process add this script to your cron tab.")
    parser.add_argument("--database",default="seqdb",help="The db prefix in the cnf file")
    parser.add_argument("--baselogdir",default="/nfs/seqscratch09/pipeline_logs/",help="Directory to log the luigi stdout and stderr, will create a timestamp specific and a sample specific directory in here, timestamp is at the level of a day")
    args = parser.parse_args()
    main(args)
