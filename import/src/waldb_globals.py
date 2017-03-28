"""
Constant variables shared across modules for the WalDB pipeline
"""
import os
import gzip
import MySQLdb
import argparse
import sys
import luigi
from operator import lt, le
from ConfigParser import RawConfigParser
import logging
from db_statements import GET_SAMPLE_DIRECTORY
from itertools import chain
from collections import OrderedDict, Counter, defaultdict
from functools import wraps
import time

cfg = RawConfigParser()
cfg.read(os.path.join(os.path.dirname(os.path.realpath(__file__)), "waldb.cfg"))
LOGGING_LEVELS = {
    "CRITICAL":logging.CRITICAL, "ERROR":logging.ERROR,
    "WARNING":logging.WARNING, "INFO":logging.INFO, "DEBUG":logging.DEBUG}
# exclude variant calls below this read depth

class DereferenceKeyAction(argparse.Action):
    """Define a class for automatically converting the key specified from
    choices to its corresponding value
    """
    def __init__(self, option_strings, dest, nargs=None, default=None,
                 choices=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs should not be specified")
        if type(choices) is not dict:
            raise TypeError("choices must be a dict")
        super(DereferenceKeyAction, self).__init__(
            option_strings, dest, choices=choices, **kwargs)
        if default:
            self.default = choices[default]

    def __call__(self, parser, namespace, values, option_string):
        setattr(namespace, self.dest, self.choices[values])

VCF_COLUMNS = ["CHROM", "POS", "rs_number", "REF", "ALT", "QUAL", "FILTER",
               "INFO", "FORMAT", "call"]
VCF_COLUMNS_DICT = dict([(column, x) for x, column in enumerate(VCF_COLUMNS)])
VALID_GTS = set(["0", "1"]) # valid values in GT field, i.e. REF/ALT
# the table format to output for calls
VARIANT_CALL_FORMAT = ("{" + "}\t{".join(
    ["sample_id", "variant_id", "block_id", "GT", "DP", "AD_REF",
     "AD_ALT", "GQ", "VQSLOD",
     "FS", "MQ", "QD", "QUAL", "ReadPosRankSum",
     "MQRankSum", "FILTER", "highest_impact", "PID_variant_id", "PGT",
     "HP_variant_id", "HP_GT", "PQ"]) + "}")
NOVEL_VARIANT_OUTPUT_FORMAT = (
    "{" + "}\t{".join(
        ["variant_id", "POS", "REF", "ALT", "rs_number", "transcript_stable_id",
         "effect_id", "HGVS_c", "HGVS_p", "polyphen_humdiv",
         "polyphen_humvar", "gene", "indel_length",
         "has_high_quality_call"]) + "}")
NOVEL_INDEL_OUTPUT_FORMAT = (
    "{" + "}\t{".join(
        ["variant_id", "POS", "REF", "ALT", "indel_length"]) + "}")
MATCHED_INDEL_OUTPUT_FORMAT = (
    "{CHROM}\t{variant_id}\t{POS}\t{REF}\t{ALT}\t{sample_id}")
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
# the regex SnpEff follows for outputting missense changes, allowing for more
# than change due to MNPs
HGVS_P_PATTERN = (
    r"^p.(?:[A-Z][a-z]{2})+(?P<codon_position>\d+)"
    r"(?P<amino_acid_changes>[A-Z][a-z]{2})+$")
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
    try:
        return MySQLdb.connect(
            read_default_file=defaults_file,
            read_default_group="client{}".format(db))
    except Exception as e:
        raise ValueError("specified database group {} is invalid in {}; error: "
                         "{}".format(db, defaults_file, e))

def get_last_insert_id(cur):
    """return the last autoincrement id for the cursor
    """
    cur.execute("SELECT LAST_INSERT_ID()")
    return cur.fetchone()[0]

def create_INFO_dict(INFO):
    """return a dict of key, value pairs in the INFO field
    """
    return dict(entry.split("=", 1) for entry in INFO.split(";") if "=" in entry)

def create_call_dict(FORMAT, call):
    """return a dict of attributes in the FORMAT field matched with the call
    """
    return dict(zip(FORMAT.split(":"), call.split(":")))

def VCF_fields_dict(line_fields, vcf_columns=VCF_COLUMNS):
    """return a dict of VCF columns for the line
    """
    return dict(zip(vcf_columns, line_fields))

def simplify_REF_ALT_alleles(REF, ALT):
    """take a potentially complex representation of a pair of alleles introduced
    into multiallelic sites and simplify to the canonical form
    http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/
    """
    # strip shared suffix
    strip = 0
    for x in xrange(1, min(len(REF), len(ALT))):
        if REF[-x] == ALT[-x]:
            strip -= 1
        else:
            break
    if strip:
        REF = REF[:strip]
        ALT = ALT[:strip]
    # strip shared prefix
    strip = 0
    for x in xrange(0, min(len(REF), len(ALT)) - 1):
        if REF[x] == ALT[x]:
            strip += 1
        else:
            break
    # return simplified REF, ALT, and position offset
    return REF[strip:], ALT[strip:], strip

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    """multiple inheritance of two argparse formatters
    """
    pass

def file_exists(arg):
    """check if the given file exists
    """
    if os.path.isfile(arg):
        return os.path.realpath(arg)
    else:
        raise argparse.ArgumentTypeError(arg + " does not exist")

def valid_numerical_argument(
    arg, arg_name, arg_type=int, min_value=0, max_value=sys.maxint,
    left_op=lt, right_op=le):
    """Confirm that the specified value is valid in the range
    (minimum_value, maximum_value] (by default)
    :param arg: the value to be tested
    :param arg_name: the name of the parameter
    :param arg_type: the type of the parameter, e.g. int or float
    :param min_value: the minimum value for the parameter, exclusive
    :param max_value: the maximum value for the parameter, inclusive
    :param left_op: the operator for testing left_op(min_value, value)
    :param right_op: the operator testing right_op(value, max_value)
    :return: arg_type(arg) if arg is valid
    """
    try:
        value = arg_type(arg)
        if left_op(min_value, value) and right_op(value, max_value):
            return value
        else:
            raise argparse.ArgumentTypeError(
                "{arg_name} ({arg}) is not in the range "
                "{left_endpoint}{min_value}, {max_value}{right_endpoint}".format(
                    arg_name=arg_name, arg=arg, min_value=min_value,
                    max_value=max_value,
                    left_endpoint="(" if left_op == lt else "[",
                    right_endpoint="]" if right_op == le else ")"))
    except TypeError:
        raise argparse.ArgumentTypeError(
            "{arg_name} ({arg}) is not a valid {arg_type}".format(
                arg_name=arg_name, arg=arg, arg_type=arg_type.__name__))

def merge_dicts(*dict_list):
    """merge an arbitrary number of dictionaries into a single dictionary and
    return it
    """
    ndicts = len(dict_list)
    if ndicts == 0:
        raise ValueError("no dicts passed in!")
    else:
        new_dict = dict_list[0].copy()
        for d in dict_list[1:]:
            new_dict.update(d)
        return new_dict

def get_data_directory(sample_name, prep_id):
    """Return the directory to the sample's data; this is stored in
    dragenQCMetrics
    """
    seqdb = get_connection("seqdb")
    try:
        seq_cur = seqdb.cursor()
        seq_cur.execute(GET_SAMPLE_DIRECTORY.format(
            prep_id=prep_id))
        row = seq_cur.fetchone()
        if row:
            path = row[0]
            if not path:
                raise OSError("Data directory for {sample_name}:{prep_id} is not "
                              "set".format(sample_name=sample_name, prep_id=prep_id))
            return os.path.join(
                path, "{sample_name}.{prep_id}".format(
                    sample_name=sample_name, prep_id=prep_id))
        else:
            raise OSError("Can't find directory for {sample_name}:{prep_id}".
                          format(sample_name=sample_name, prep_id=prep_id))
    finally:
        if seqdb.open:
            seqdb.close()

def strip_prefix(string, prefix):
    """Strip the given prefix if it's present
    """
    return string[len(prefix):] if string.startswith(prefix) else string

def strip_suffix(string, suffix):
    """Strip the given suffix if it's present
    """
    return string[:-len(suffix)] if string.endswith(suffix) else string

class MultipleFilesTarget(luigi.Target):
    """Target which checks for existence of all files passed in
    """
    def __init__(self, targets):
        self.targets = list(chain(*targets))

    def exists(self):
        for fn in self.targets:
            if not os.path.isfile(fn):
                return False
        return True

    def get_targets(self):
        return self.targets

# keep a count for each key and maintain the order in which keys were added
class OrderedCounter(Counter, OrderedDict):
    pass

# initialize any new key to the default value and maintain the order in which
# keys were added
class OrderedDefaultDict(OrderedDict, defaultdict):
    def __init__(self, default_factory=None, *args, **kwargs):
        super(OrderedDefaultDict, self).__init__(*args, **kwargs)
        self.default_factory = default_factory

# wrapper to time a function
def timer(fh=sys.stdout):
    def timer_wrapper(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            current_time = time.time()
            call = ("{func}({args}{separator}{kwargs})".format(
                func=func.func_name, args=", ".join([str(arg) for arg in args]),
                separator=", " if args and kwargs else "",
                kwargs=", ".join(["{}={}".format(arg, value) for arg, value in
                                  kwargs.iteritems()])))
            fh.write("Calling {call} @{ctime}\n".format(call=call, ctime=time.ctime()))
            result = func(*args, **kwargs)
            fh.write("Finished {call} @{ctime}; elapsed time: {elapsed_time} seconds\n".
                  format(call=call, ctime=time.ctime(),
                         elapsed_time=time.time() - current_time))
            return result
        return wrapper
    return timer_wrapper