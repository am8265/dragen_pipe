#!/nfs/goldstein/software/python2.7.7/bin/python

import MySQLdb
import os
import argparse
import subprocess
import glob 
import sys
import time

global  luigi_tasks = ['RunCvgMetrics','Binning','AlignmentMetrics','DuplicateMetrics',
                   'VariantCallingMetrics','CheckUpdates','UpdateDatabase','EthnicityPipeline','RunRelatednessCheck']

global jobs

def check_basedir(basedir):
    try :
        basedir=os.path.abspath(basedir) ## Strip of terminal '/' for uniformity 
        return basedir
    except:
        print "Malformed base directory !!\nPlease check again...\n\n"
        return 'malformed'

def retrieve_localsamples(basedir):
    """ Function to retrieve paths for all available samples to be processed
    """

    for f in glob.glob(basedir+'/*'):
        yield f


def is_exists_file(search_dir,sample,ftype):
    """
    Check whether a gvcf/bam with the sample stem is present in 
    the directory, if ftype is strict , check for both
    """

    b = os.path.isfile(search_dir+'/{0}.realn.recal.bam'.format(sample))
    v = os.path.isfile(search_dir+'/{0}.g.vcf.gz'.format(sample))
    if ftype == 'bam':
        return b
    if ftype == 'gvcf':
        return v
    if ftype == 'strict':
        return (b and v)
    
def map_tasks_to_luigi():
    """ 
    Map user entered tasks to luigi tasks
    return type : list 
    """
    
def check_running():
    """
    """

    for p in range(len(jobs):0:-1): ## Iterating from the end of the list 
        if jobs[p].poll() is not None: # A none value indicates that the process has not yet terminated
            del jobs[p] ## Remove if done 
            
def wrapper_local(args):
    """ Running a local version
    """
    
    print "Running Local Mode...\n\n"

    ## Ensure the basedirectory is proper
    if check_basedir(args.basedir) == 'malformed':
        sys.exit()
    else:
        args.basedir = check_basedir(args.basedir)

    ## Retrieve the all the sample directories
    sample_dirs = retrieve_localsamples(args.basedir)
    
    for s in sample_dirs:
        samplename=s.split('/')[-1]
        if is_exists_files(s,samplename,'strict'): ## Relaxed criteria for now
            for t in args.task:          
                proc = run_cvgbin(t,samplename,basedir)
                jobs.append(proc)
                if len(jobs) < 10:
                    jobs.append(proc)
                else:
                    time.sleep(100) ## Tinker this value 
                    check_running()

                    
def wrapper_pipeline(args):
    """ The pipeline version
    """
    print "Running in Pipeline Mode...\n\n"
    
        
def main(args):
    
    args.func(args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser('A Wrapper for the qc modules',description='Can be run either in the pipeline or local mode, this is a work in progress')
    parser.add_argument('--all','-all',help='Run all tasks',action='store_true',
                        default=False)
    parser.add_argument('--task','-task',help='Specify task(s)',choices=luigi_tasks,required=True)
    subparsers = parser.add_subparsers(help='sub-command')
    parser_pipeline = subparsers.add_parser('Pipeline',help='Mode of operation, i.e. Pipeline or Local')
    parser_pipeline.add_argument('-cnf','--cnf',help='MySQL config file')
    parser_pipeline.add_argument('-testdb','--testdb',help='Run on testdb',action='store_true',
                        default=False)

    parser_pipeline.add_argument('--update','-update',help='Whether to update to database',action='store_true',default='False')

    parser_pipeline.set_defaults(func=wrapper_pipeline)

    parser_local = subparsers.add_parser('Local',help='Mode of operation, i.e. Pipeline or Local')
    
    parser_local.add_argument('--basedir','-basedir',help='Base directory from which to look for samples\nE.g. /nfs/seqscratch11/Exomes/',required=True)
    
    parser_local.set_defaults(func=wrapper_local)
    args = parser.parse_args()
    main(args)
