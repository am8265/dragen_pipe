import MySQLdb
import os
from ConfigParser import ConfigParser
import subprocess
import sys
import luigi
from luigi.contrib.sge import SGEJobTask
from datetime import datetime

def get_familyid(sample_name):
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

##### pointless...
class MyExtTask(luigi.ExternalTask):
    file_loc = luigi.Parameter()
    def output(self):
        return luigi.LocalTarget(self.file_loc)

def get_connection(config_file="/nfs/goldstein/software/dragen_pipe/master/automated_ethnicity_relatedness/database.cfg",db="seqdb"):
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

# wtf?!? no reason for dps 31 - i.e. it's just off of alignseqfileloc?!
# """ SELECT Q.CHGVID,Q.pseudo_prepid,Q.SeqType,Q.AlignSeqFileLoc FROM"""
    query = ("SELECT sample_name, d.pseudo_prepid, sample_type from dragen_sample_metadata d \
      JOIN dragen_ped_status p on d.pseudo_prepid=p.pseudo_prepid \
      where p.create_ped = 1 and p.append_ped = 0 order by pseudo_prepid desc" )

    print("using '{}'".format(query))
    db = get_connection(db="seqdb")
    cur = db.cursor()
    cur.execute(query)
    result = cur.fetchall()
    return result

# f'ing pointless - CreatePed is dependency of PredictAndUpdate and sets create_ped while PredictAndUpdate does predict and AppendMasterPed does append_ped
def update_ped_status(pseudo_prepid,field):
    """ Update status 

    pseudo_prepid : str 
    field : str
    """
    update_statement = ("UPDATE dragen_ped_status SET {0} = 1 WHERE pseudo_prepid = {1}""".format(
      field,pseudo_prepid ) )    
    db = get_connection(db="seqdb")
    try:
        cur = db.cursor()
        cur.execute(update_statement)
        db.commit()
    finally:
        db.close()

def run_shellcmd(cmd):
    proc = subprocess.Popen(cmd,shell=True)
    proc.wait()
    if proc.returncode: ## Non zero return code
        raise Exception( subprocess.CalledProcessError( proc.returncode, cmd ) )
        
