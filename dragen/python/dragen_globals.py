"""
Constant variables shared across modules for the DRAGEN pipeline
"""

import MySQLdb
import os
from ConfigParser import ConfigParser

cfg = ConfigParser()
cfg.read(os.path.join(os.path.dirname(os.path.realpath(__file__)), "luigi.cfg"))

def get_connection(db):
    """return a connection to the database specified
    """

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

def get_cfg():
    return cfg
