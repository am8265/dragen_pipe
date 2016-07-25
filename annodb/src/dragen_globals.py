"""
Constant variables shared across modules for the DRAGEN pipeline
"""
import os
import gzip
import MySQLdb
from ConfigParser import ConfigParser

cfg = ConfigParser()
cfg.read(os.path.join(os.path.dirname(os.path.realpath(__file__)), "anno.cnf"))

CNF = "/nfs/goldstein/software/dragen/dragen.cnf" # defaults file for pipeline
VCF_COLUMNS = ["CHROM", "POS", "rs_number", "REF", "ALT", "QUAL", "FILTER",
               "INFO", "FORMAT", "call"]
VCF_COLUMNS_DICT = dict([(column, x) for x, column in enumerate(VCF_COLUMNS)])
BLOCK_SIZE = 10000 # the bases in a block of variant calls (for indexing)
VALID_GTS = set(["0", "1"]) # valid values in GT field, i.e. REF/ALT
# the table format to output for calls
VARIANT_CALL_FORMAT = ("{" + "}\t{".join(
    ["sample_id", "variant_id", "block_id", "GT", "DP_pileup", "DP", "AD_REF",
     "AD_ALT", "GQ", "FS", "MQ", "QD", "QUAL", "ReadPosRankSum", "MQRankSum",
     "PASS"]) + "}")
POLYPHEN_ATTRIB_ID = {"humvar":268, "humdiv":269}
# PolyPhen scores are packed by sorted one-letter codes
AMINO_ACIDS = dict(
    [[aa, aa_idx] for aa_idx, aa in enumerate((
        "Ala", #A
        "Cys", #C
        "Asp", #D
        "Glu", #E
        "Phe", #F
        "Gly", #G
        "His", #H
        "Ile", #I
        "Lys", #K
        "Leu", #L
        "Met", #M
        "Asn", #N
        "Pro", #P
        "Gln", #Q
        "Arg", #R
        "Ser", #S
        "Thr", #T
        "Val", #V
        "Trp", #W
        "Tyr" #Y
    ))])
# the regex SnpEff follows for outputting missense changes
HGVS_P_PATTERN = (r"p\.[A-Z][a-z]{2}(?P<codon_position>\d+)"
                  "(?P<amino_acid_change>[A-Z][a-z]{2})")
POLYPHEN_PROB_BITMASK = 2 ** 10 - 1

def get_fh(fn, mode="r"):
    """return a file handle to the file, gzipped optional
    """
    try:
        with open(fn) as fh:
            magic_number = fh.read(2)
    except:
        if os.path.splitext(fn)[-1] == ".gz":
            magic_number = "\x1f\x8b"
        else:
            magic_number = None
    if magic_number == "\x1f\x8b":
        return gzip.open(fn, mode)
    else:
        return open(fn, mode)

def get_cfg():
    return cfg

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

def get_last_insert_id(cur):
    """return the last autoincrement id for the cursor
    """
    cur.execute("SELECT LAST_INSERT_ID()")
    return cur.fetchone()[0]

def create_INFO_dict(INFO, call):
    """return a dict of attributes in the INFO field matched with the call
    """
    return dict(zip(INFO.split(":"), call.split(":")))
