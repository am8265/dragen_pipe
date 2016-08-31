"""
Constant variables shared across modules for the DRAGEN pipeline
"""

import MySQLdb
import os
from ConfigParser import ConfigParser
import subprocess

cfg = ConfigParser()

### Get luigi config file from the environment variable 
proc = subprocess.Popen("echo $LUIGI_CONFIG_PATH",shell=True,stdout=subprocess.PIPE)
proc.wait()
if proc.returncode: ## Non zero return code
    raise subprocess.CalledProcessError(proc.returncode,"echo $LUIGI_CONFIG_PATH\nCould not find luigi config file, check your path!\n")

config_file = proc.stdout.read().strip('\n')
cfg.read(config_file)


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
