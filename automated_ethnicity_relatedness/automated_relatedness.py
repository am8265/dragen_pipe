import MySQLdb
import os
from ConfigParser import ConfigParser
import subprocess
import sys
import luigi
from luigi.contrib.sge import SGEJobTask
from datetime import datetime

def get_familyid(sample_name):
    """
    """
    
    query = "SELECT FamilyID FROM SampleT WHERE CHGVID = '{0}'".format(sample_name)
    db = get_connection(db="seqdb")
    try:
        cur = db.cursor()
        cur.execute(query)
        result = cur.fetchall()
        familyid = result[0][0]
        if familyid == "N/A":
            return sample_name
        else:
            return familyid
    finally:
        db.close()

class MyExtTask(luigi.ExternalTask):
    """Check for files
    """

    file_loc = luigi.Parameter()
    def output(self):
        return luigi.LocalTarget(self.file_loc)

def get_connection(config_file="/nfs/seqscratch09/rp2801/ethnicity_predictions/scripts/automated_ethnicity/database.cfg",db="seqdb"):
    """ Return a connection to seqdb
    """
    
    try:
        cfg = ConfigParser()
        cfg.read(config_file)
    except Exception as e:
        raise Exception("Failed to read config file : %s"%e)

    defaults_file = cfg.get("db","cnf")
    if db == "seqdb":
        defaults_group = "client"+cfg.get("db","seqdb_group")
    else:
        raise Exception("Use db = seqdb to connect to sequencedb !")
    
    return MySQLdb.connect(
        read_default_file=defaults_file,
        read_default_group=defaults_group)


def get_samples_append():
    """
    Get samples which have been archived,
    i.e. pipeline_step_id = 31 has finished = 1
    which need to be appended to the masterped 

    Returns a nested tuple :
    ((chgvid,pseudo_prepid,seqtype,alignseqfileloc),(...),(...))
    """

    query = (
        """ SELECT Q.CHGVID,Q.pseudo_prepid,Q.SeqType,Q.AlignSeqFileLoc FROM"""
        """ dragen_qc_metrics as Q INNER JOIN dragen_pipeline_step as D"""
        """ INNER JOIN dragen_ped_status as P"""
        """ ON Q.pseudo_prepid = D.pseudo_prepid AND"""
        """ P.pseudo_prepid = D.pseudo_prepid AND"""
        """ D.pipeline_step_id = 31 AND D.finished = 1 AND"""
        """ P.append_ped = 0"""
        )
    db = get_connection(db="seqdb")
    cur = db.cursor()
    cur.execute(query)
    result = cur.fetchall()
    return result

def update_ped_status(pseudo_prepid,field):
    """ Update status 

    pseudo_prepid : str 
    field : str
    """
    update_statement = (""" UPDATE dragen_ped_status SET {0} = 1"""
                        """ WHERE pseudo_prepid = {1}""".format(
                            field,pseudo_prepid)
                        )    
    db = get_connection(db="seqdb")
    try:
        cur = db.cursor()
        cur.execute(update_statement)
        db.commit()
    finally:
        db.close()

def run_shellcmd(cmd):
    """Use subprocess to run commands
    """

    proc = subprocess.Popen(cmd,shell=True)
    proc.wait()
    if proc.returncode: ## Non zero return code
        raise Exception(subprocess.CalledProcessError(proc.returncode,
                                                                  cmd))
        
