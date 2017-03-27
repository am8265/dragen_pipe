#!/nfs/goldstein/software/python2.7.7/bin/python2.7

import argparse
import glob
import MySQLdb
import os
import shlex
import string
import sys
import subprocess
import signal
import multiprocessing as mp 
from time import sleep,time
from threading import Thread
import random
from dragen_globals import get_connection
from db_statements import * 
from datetime import datetime


##############################################################################
###                   Automates submission of the gatk                     ###
###                   luigi pipeline to the cluster                        ###                    
##############################################################################

## Define globals for storing tasks
## Tasks are keyed by their name {sample_name}.{pseudo_prepid}.{sample_type}
RUNNING = {}
EXECUTED = {}
ERRORS = {}
email = None
num_processes = 0

def execute_query(database,query,fetch):
    """ Execute db fetch statements

    database : The database to use for the query
    query : The query to be performed
    fetch : Bool; Indicates whether to fetch the results incase of a SELECT statement
    
    return : Nothing
    """
    try:
        db = get_connection(database)
        cur = db.cursor()
        cur.execute(query)
        if fetch == True: ## For select statements 
            result = cur.fetchall()
            return result
        else: ## For insert and delete statements
            db.commit()
    except MySQLdb.Error, e:
        clean_tasks(force_exit=True)
        raise Exception("Error %d IN CONNECTION: %s" %(e.args[0],e.args[1]))
    finally:
        db.close()
            
def get_samples():
    """ Get samples from the gatk_queue table
    """
    query = GET_SAMPLES_TO_RUN
    result = execute_query("seqdb",query,fetch=True)
    return result

def create_logdir(baselogdir,sample_name,pseudo_prepid):
    """ Create luigi stdout and stderr logging directory for a sample
    
    baselogdir : str ; Base logging directory
    sample_name : str ; The sample name
    pseudo_prepid : int ; The pseudo_prepid

    return : string ; Final path to the sample log directory 
    """
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d")
    sample_log_dir = os.path.join(timestamp+'_'+'logs',sample_name+'_'+str(pseudo_prepid))
    final_log_dir = os.path.join(baselogdir,sample_log_dir)
    if not os.path.exists(final_log_dir):
        os.makedirs(final_log_dir)
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
         result = execute_query("seqdb",query,fetch=True)
         return result[0][0]
    else:
        query = ("SELECT region_file_lsrc "
                 "FROM captureKit WHERE prepT_name='{0}' "
                 "AND chr = 'all'").format(capture_kit)
        result = execute_query("seqdb",query,fetch=True)
        return result[0][0]

def clean_tasks(force_exit = False,test=False):
    """ Clean up luigi pipeline commands
    """
    for name, task_details in RUNNING.items():
        if not task_details[0].is_alive():
            if task_details[0].exitcode:## Non zero return code, the task failed
                ERRORS[name] = task_details[0].exitcode
            else: ## Success
                EXECUTED[name] = 'success!'
                RUNNING.pop(name)
        else: ## Still running , check if the running_time has exceeded
            cur_time = datetime.now()
            start_time = task_details[1]
            diff = cur_time - start_time
            time_elapsed = diff.total_seconds()
            sample_type = name.split('.')[-1]
            if force_exit == True:
                print "Exiting wrapper; Killing sample %s ....\n"%name
                kill_process(task_details[0])
                kill_qstat_jobs(name, test = test)
            else:
                if sample_type.upper() == 'GENOME':
                    if time_elapsed > 518400: ## Exceeded 6 days
                        print "Genome sample %s exceeded 6 days Killing process ...\n"%name
                        kill_process(task_details[0])
                        kill_qstat_jobs(name,test=test)
                        RUNNING.pop(name)
                else:
                    if time_elapsed > 54000: ## Exceeded 15 hrs
                        print "Exome sample %s exceeded 15 hours Killing process ...\n"%name
                        kill_process(task_details[0])
                        kill_qstat_jobs(name,test=test)
                        RUNNING.pop(name)
   
def kill_process(process):
    """
    process : a multiprocessing process object
    """
    pid = process.pid
    os.kill(pid,signal.SIGKILL)

def kill_qstat_jobs(task_name,test=False):
    """
    task_name : str; e.g. sample_name.pseudo_prepid.sample_type
    test : bool ; whether to run in test mode 
    """
    if test:
        return True
    
    ## Since qstat jobs are defined as sample_name.pseudo_prepid in gatk_pipeline
    job_name = '.'.join(task_name.split('.')[0:2])
    shell_cmd = ("""qstat -r -ne | grep -B 1 %s|"""
                 """cut -f1 -d' '|xargs qdel -j {}"""%(job_name))
    proc = subprocess.Popen(shell_cmd,shell=True)
    proc.wait()

def test_luigi(sample_info,baselogdir):
    """ This is a test target function for testing purposes
    """
    pseudo_prepid,sample_name,sample_type,capture_bed,scratch = sample_info
    final_log_dir = create_logdir(baselogdir,sample_name,pseudo_prepid)
    out_file = os.path.join(final_log_dir,sample_name+'_'+pseudo_prepid+'.stdout')
    error_file = os.path.join(final_log_dir,sample_name+'_'+pseudo_prepid+'.stderr')

    try:
        flag = 0
        key = sample_name+'.'+pseudo_prepid+'.'+sample_type
        cmd = (
            """ echo luigi --module gatk_pipe ArchiveSample --sample-name {0} --pseudo-prepid {1} --capture-kit-bed  {2} --sample-type {3} --scratch {4} --poll-time 120 --workers 2 --worker-wait-interval 180 --scheduler-remove-delay 86400 >> {5} 2>>{6} & sleep 100""".format(sample_name,pseudo_prepid,capture_bed,sample_type,scratch,out_file,error_file))
        proc = subprocess.Popen(cmd, shell = True)
        proc.wait()
        exit(proc.returncode)
    except KeyboardInterrupt:
        flag = 1
    finally:
        if flag == 1:
            print "Keyboard interrupt, cleaning up process ... {0}.{1} \n".format(sample_name,pseudo_prepid)
    
def run_luigi(sample_info,baselogdir):
    """ Execute the luigi command with the appropriate parameters as a system call     

    sample_info : tuple ; (pseudo_prepid,sample_name,...)
    baselogdir : str; base directory to store the logs ; e.g. /nfs/seqscratch09/luigi_logs_test/
    """   
    pseudo_prepid,sample_name,sample_type,capture_bed,scratch = sample_info
    final_log_dir = create_logdir(baselogdir,sample_name,pseudo_prepid)
    out_file = os.path.join(final_log_dir,sample_name+'_'+pseudo_prepid+'.stdout')
    error_file = os.path.join(final_log_dir,sample_name+'_'+pseudo_prepid+'.stderr')
    
    try:
        flag = 0
        key = sample_name+'.'+pseudo_prepid+'.'+sample_type
        #cmd = (""" luigi --module gatk_pipe ArchiveSample --sample-name {0} --pseudo-prepid {1} --capture-kit-bed  {2} --sample-type {3} --scratch {4} --poll-time 120 --workers 2 --worker-wait-interval 180 --scheduler-remove-delay 86400 >> {5} 2>>{6}""".format(sample_name,pseudo_prepid,capture_bed,sample_type,scratch,out_file,error_file))
        cmd = (""" luigi --module gatk_pipe ArchiveSample --sample-name {0} --pseudo-prepid {1} --capture-kit-bed  {2} --sample-type {3} --scratch {4} --poll-time 120 --workers 2 --worker-wait-interval 180 --scheduler-remove-delay 86400""".format(sample_name,pseudo_prepid,capture_bed,sample_type,scratch))
        with open(out_file,'w') as OUT, open(error_file,'w') as ERR:
            proc = subprocess.Popen(shlex.split(cmd),stdout=OUT,stderr=ERR)
            proc.wait()
            exit(proc.returncode)            
    except KeyboardInterrupt:
        flag = 1
    finally:
        if flag == 1:
            print "Keyboard interrupt, cleaning up process ... {0}.{1} \n".format(sample_name,pseudo_prepid)
        
def start_pipeline(sample_info,baselogdir,test=False):
    """ Start a luigi gatk pipeline     
    sample_info : tuple ; parameters for running a luigi job
    baselogdir : str ; the base directory to store the logs
    test : bool ; whether to run in test mode
    """
    pseudo_prepid,sample_name,sample_type,capture_bed,scratch = sample_info
    name = sample_name+'.'+pseudo_prepid+'.'+sample_type
    if name not in RUNNING:
        if test:
            p = mp.Process(target=test_luigi,args=[sample_info,baselogdir],name=name)
        else:
            p = mp.Process(target=run_luigi,args=[sample_info,baselogdir],name=name)
        RUNNING[name] = (p,datetime.now())
        p.start()
        
   
def main(args):
    """ The main function
    """    
    ## Unpack parameters 
    ## Dont want to pass email around to clean_tasks,keeping it global should be pretty safe
    global email
    email = args.email
    baselogdir = args.baselogdir
    max_processes = args.max_processes
    wait_time = args.wait_time
    scratch = "/nfs/%s/ALIGNMENT/BUILD37/DRAGEN"%args.scratch
    test = args.test

    try:
        while True:
            for pseudo_prepid,sample_name,sample_type,capture_kit,priority in get_samples():
                task_key = sample_name+'.'+str(pseudo_prepid)+'.'+sample_type
                if task_key not in RUNNING and task_key not in ERRORS:
                    capture_bed = get_capturekit_info(sample_type,capture_kit)
                    sample_tup = (str(pseudo_prepid),sample_name,sample_type,capture_bed,scratch)
                    wild_card_bam = glob.glob(scratch+'/%s/%s.%s/%s.%s*.bam'%(sample_type.upper(),sample_name,pseudo_prepid,sample_name,pseudo_prepid))
                    if len(RUNNING.keys()) <= max_processes:
                        if len(wild_card_bam)>0:
                            start_pipeline(sample_tup,baselogdir,test=test)
                        else:
                            print "Will not run {0}.{1} bam not found !".format(sample_name,pseudo_prepid)
                    else:
                        while len(RUNNING.keys()) >= max_processes: ## Hold till processes are freed up
                            print "Tasks Executed : ",EXECUTED
                            print "Tasks Running : ",RUNNING
                            print "Failed Tasks :",ERRORS
                            print "Task Queued : ",task_key                            
                            sleep(wait_time)
                            clean_tasks(force_exit = False,test=test)            
    except KeyboardInterrupt:
        print "Keyboard interrupt in main"
        for key in RUNNING:
            kill_process(RUNNING[key][0])
            kill_qstat_jobs(key)
    finally:
        print "Exiting Main....\n"
        sys.exit()                          
        
                               
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Luigi wrapper script',description="This script is a wrapper for running the new dragen variant calling and qc pipeline.")
    parser.add_argument("--baselogdir",default="/nfs/seqscratch09/pipeline_logs/",help="Directory to log the luigi stdout and stderr, will create a timestamp and sample specific directory in here, timestamp is at the level of a day",required = False)
    parser.add_argument("--max-processes",dest='max_processes',default=300,help="The maximum number of processes to run",required=False,type=int)
    parser.add_argument("--wait-time",dest='wait_time',default=3600,help="The number of seconds to wait before querying the database for samples again",required=False,type=int)
    parser.add_argument("--scratch",default="seqscratch_ssd",help="The scratch directory to operate on",required=False)
    parser.add_argument("--email",default=None,help="Email to send notifications to",required=False)
    parser.add_argument("--test",action='store_true',help="Will simply echo the luigi commands instead of running them")
    args = parser.parse_args()
    main(args)
