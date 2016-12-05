"""
Constant variables shared across modules for the DRAGEN pipeline
"""

import MySQLdb
import os
from ConfigParser import ConfigParser
import subprocess


def get_connection(db):
    """return a connection to the database specified
    """
    cfg = ConfigParser()

    ## Look in the current directory for the database config file 
    config_file = r"/home/rp2801/git/dragen_pipe/dragen/python/database.cfg"
    try:
        cfg.read(config_file)
    except Exception as e :
        print str(e)
        
    cfg.read(config_file)

    defaults_file = cfg.get("db", "cnf")
    if db == "dragen":
        return MySQLdb.connect(
            read_default_file=defaults_file,
            read_default_group="client" + cfg.get("db", "dragen_group"))
    elif db == "seqdb":
        return MySQLdb.connect(
            read_default_file=defaults_file,
            read_default_group="client" + cfg.get("db", "seqdb_group"))
    else:
        raise ValueError("specified database group is invalid")
