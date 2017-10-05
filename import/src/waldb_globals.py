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
from ConfigParser import ConfigParser, RawConfigParser
import logging
from db_statements import (
    GET_SAMPLE_DIRECTORY, GET_TIMES_STEP_RUN, BEGIN_STEP,
    FINISH_STEP, FAIL_STEP, GET_PIPELINE_STEP_ID, GET_STEP_STATUS,
    GET_SAMPLE_METADATA, GET_CAPTURE_KIT_BED, INSERT_PIPELINE_STEP)
from itertools import chain
from collections import OrderedDict, Counter, defaultdict
from functools import wraps
import time
import subprocess
from luigi.contrib.sge import SGEJobTask
from shlex import split as sxsplit
from pprint import pprint

cfg = RawConfigParser()
cfg.read(os.path.join(os.path.dirname(os.path.realpath(__file__)), "waldb.cfg"))
LOGGING_LEVELS = {
    "CRITICAL":logging.CRITICAL, "ERROR":logging.ERROR,
    "WARNING":logging.WARNING, "INFO":logging.INFO, "DEBUG":logging.DEBUG}
CHROMs = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
          "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"]

def get_pipeline_version():
    version = subprocess.check_output(
        ["/nfs/goldstein/software/git-2.5.0/bin/git", "describe", "--tags"]).strip()
    if version:
        return version
    else:
        raise ValueError("Could not get the version # of the pipeline; "
                         "maybe run it from a directory in the repo?")

class SQLTarget(luigi.Target):
    """ A luigi target class describing verification of the entries in the database
    """
    def __init__(self, pseudo_prepid, pipeline_step_id):
        self.pseudo_prepid = pseudo_prepid
        self.pipeline_step_id = pipeline_step_id

    def exists(self):
        db = get_connection("seqdb")
        try:
            cur = db.cursor()
            cur.execute(GET_STEP_STATUS.format(
                prep_id=self.pseudo_prepid,
                pipeline_step_id=self.pipeline_step_id))
            row = cur.fetchone()
            if row:
                return row[0] == "completed"
            else:
                cur.execute(INSERT_PIPELINE_STEP.format(
                    prep_id=self.pseudo_prepid,
                    pipeline_step_id=self.pipeline_step_id))
                db.commit()
                return False
        finally:
            if db.open:
                db.close()

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

def run_command(command, directory, task_name, print_command=True):
    """Run the specified command in a subprocess and log the output
    """
    if print_command:
        print(command)
    base_name = os.path.join(directory, task_name)
    with open(base_name + ".out", "w") as out_fh, \
            open(base_name + ".err", "w") as err_fh:
        out_fh.write(command + "\n")
        out_fh.flush()
        p = subprocess.Popen(sxsplit(command), stdout=out_fh, stderr=err_fh)
        p.wait()
    if p.returncode:
        raise subprocess.CalledProcessError(p.returncode, command)

def file_handle(f, mode="w"):
    """Return a file handle to f; if it's a file handle already, return it,
    otherwise open it
    """
    if type(f) is str:
        return open(f, mode)
    if type(f) is file:
        return f
    if type(f) is int:
        return f
    raise TypeError("Wrong argument passed in for opening")

def close_file_handles(file_handles):
    """Close the file handles as appropriate
    """
    for file_handle in file_handles:
        if (type(file_handle) is file and not file_handle.closed
            and file_handle not in (sys.stdout, sys.stderr)):
            file_handle.close()

class PipelineTask(SGEJobTask):
    """Abstract Task class for automatically updating pipeline step statuses
    """
    pseudo_prepid = luigi.IntParameter(
        description="The pseudo_prepid for this sample; used for "
        "obtaining/updating statuses")
    print_init = luigi.BoolParameter(default=False)

    def __init__(self, *args, **kwargs):
        super(PipelineTask, self).__init__(*args, **kwargs)
        self.pipeline_step_id = self._get_pipeline_step_id()
        # these will be passed onto subprocess, overwrite in pre_shell_commands
        # method if needed
        self.shell_options = {"record_commands_fn":None, "stdout":os.devnull,
                              "stderr":None, "shell":False}
        self.commands = []
        self.files = []
        self.directories = []
        if self.print_init:
            print("Initializing {}".format(self.__class__.__name__))

    def _get_pipeline_step_id(self):
        """Set the pipeline_step_id for any class inheriting from this
        """
        seqdb = get_connection("seqdb")
        try:
            seq_cur = seqdb.cursor()
            seq_cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            row = seq_cur.fetchone()
            if row:
                return row[0]
            else:
                raise ValueError("Could not find pipeline step of name {}!".
                                 format(self.__class__.__name__))
        finally:
            if seqdb.open:
                seqdb.close()

    def _set_pipeline_version(self):
        """Set the version of the pipeline being executed by checking the git
        respository's version tag
        """
        self.version = get_pipeline_version()

    def run(self):
        """First create/update as needed the record in the pipeline step table,
        set its start time, and increment its times_ran counter
        """
        self._set_pipeline_version()
        self._update_step_start_time()
        try:
            super(PipelineTask, self).run()
            self._update_step_success()
            self._run_post_success()
        except:
            self._update_step_failure()
            raise

    def _update_step_start_time(self):
        seqdb = get_connection("seqdb")
        try:
            seq_cur = seqdb.cursor()
            seq_cur.execute(GET_TIMES_STEP_RUN.format(
                prep_id=self.pseudo_prepid,
                pipeline_step_id=self.pipeline_step_id))
            row = seq_cur.fetchone()
            times_ran = row[0] if row else 0
            seq_cur.execute(BEGIN_STEP.format(
                pseudo_prepid=self.pseudo_prepid,
                pipeline_step_id=self.pipeline_step_id,
                version=self.version, times_ran=times_ran + 1))
            seqdb.commit()
        finally:
            if seqdb.open:
                seqdb.close()

    def _update_step_success(self):
        seqdb = get_connection("seqdb")
        try:
            seq_cur = seqdb.cursor()
            seq_cur.execute(FINISH_STEP.format(
                pseudo_prepid=self.pseudo_prepid,
                pipeline_step_id=self.pipeline_step_id))
            seqdb.commit()
        finally:
            if seqdb.open:
                seqdb.close()

    def _run_post_success(self):
        """Add anything here that should be performed only after successful
        update of the database with the success of the task, e.g. deleting files
        """
        pass

    def _update_step_failure(self):
        seqdb = get_connection("seqdb")
        try:
            seq_cur = seqdb.cursor()
            seq_cur.execute(FAIL_STEP.format(
                pseudo_prepid=self.pseudo_prepid,
                pipeline_step_id=self.pipeline_step_id))
            seqdb.commit()
        finally:
            if seqdb.open:
                seqdb.close()
                 
    def work(self):
        """Always execute the following three steps in order:
            1. pre-shell commands
            2. shell commands
            3. post-shell commands
            4. updating relevant files/directories to group bioinfo and
            permissions to 664 and 775 for files and directories, respectively
        Don't overwrite this method, instead overwrite pre- and post-shell
        command functions as needed
        """
        os.umask(0002)
        self.pre_shell_commands()
        self.run_shell_commands()
        self.post_shell_commands()
        self.update_permissions()

    def pre_shell_commands(self):
        """Overwrite to execute anything code that should run prior to running
        shell comamnds
        """
        pass

    def post_shell_commands(self):
        """Overwrite to execute any code that should run after running shell
        commands
        """
        pass

    def run_shell_commands(self):
        """Run an arbitrary list of commands and raise a CalledProcessError in
        the event that any returns a non-zero error code
        """
        with open(self.shell_options.pop("record_commands_fn"), "w") as out:
            out.write("#Pipeline version #:{version}\n".format(version=self.version))
            out.write("\n".join(self.commands) + "\n")
        if self.commands:
            for command in self.commands:
                if not command:
                    continue
                fhs = {"stdout":None, "stderr":None}
                try:
                    if not self.shell_options["shell"]:
                        command = sxsplit(command)
                    if not self.shell_options["stderr"]:
                        fhs["stderr"] = file_handle(sys.stdout)
                    fh_dict = {}
                    for fh in ("stdout", "stderr"):
                        f = self.shell_options.pop(fh)
                        # store these in a temporary dict and add back later, as
                        # subprocess will raise an exception if these are passed
                        # in
                        fh_dict[fh] = f
                        if f:
                            fhs[fh] = file_handle(f)
                    p = subprocess.Popen(
                        command, stdout=fhs["stdout"], stderr=fhs["stderr"],
                        **self.shell_options)
                    p.wait()
                    if p.returncode:
                        raise subprocess.CalledProcessError(
                            p.returncode, command)
                    self.shell_options.update(fh_dict)
                finally:
                    close_file_handles([fhs["stdout"], fhs["stderr"]])

    def update_permissions(self):
        for fn in self.files:
            os.chmod(fn, 0664)
        for d in self.directories:
            os.chmod(d, 0775)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

class DragenPipelineTask(PipelineTask):
    """Add the parameters that are de facto in every gatk_pipe task to reduce
    repetitiveness
    """
    sample_name = luigi.Parameter(
        default=None, description="the sample identifier")
    capture_kit = luigi.Parameter(
        default=None, description="the capture kit used")
    capture_kit_bed = luigi.InputFileParameter(
        default=None,
        description="the location of the BED file containing the capture kit regions")
    sample_type = luigi.ChoiceParameter(
        default=None, choices=["EXOME", "GENOME"],
        description="The type of sequencing performed for this sample")
    scratch = luigi.Parameter(
        default=None,
        description="The scratch space to use for the pipeline")

    def __init__(self, *args, **kwargs):
        super(DragenPipelineTask, self).__init__(*args, **kwargs)
        if (not self.sample_name or not self.capture_kit_bed
            or not self.sample_type or not self.scratch):
            seqdb = get_connection("seqdb")
            try:
                seq_cur = seqdb.cursor()
                seq_cur.execute(
                    GET_SAMPLE_METADATA.format(prep_id=self.pseudo_prepid))
                row = seq_cur.fetchone()
                if row:
                    sample_name, sample_type, capture_kit, scratch, _ = row
                    if not self.sample_name:
                        kwargs["sample_name"] = sample_name
                    if not self.sample_type:
                        kwargs["sample_type"] = sample_type.upper()
                    if not self.scratch:
                        kwargs["scratch"] = scratch
                    if not self.capture_kit:
                        kwargs["capture_kit"] = capture_kit
                    if not self.capture_kit_bed:
                        seq_cur.execute(GET_CAPTURE_KIT_BED.format(
                            capture_kit="Roche"
                            if kwargs["sample_type"] == "GENOME" else
                            capture_kit))
                        kwargs["capture_kit_bed"] = seq_cur.fetchone()[0]
                else:
                    raise ValueError("Couldn't find sample metadata for {}!".
                                     format(self.pseudo_prepid))
                super(DragenPipelineTask, self).__init__(*args, **kwargs)
            finally:
                if seqdb.open:
                    seqdb.close()

class GATKPipelineTask(DragenPipelineTask):
    """Add implicit parameters that will always be found in specific paths for
    data at various stages of the pipeline
    """
    def __init__(self, *args, **kwargs):
        super(GATKPipelineTask, self).__init__(*args, **kwargs)
        name_prep = "{}.{}".format(self.sample_name, self.pseudo_prepid)
        self.name_prep = name_prep
        self.scratch_dir = os.path.join(
            "/nfs", self.scratch, "ALIGNMENT", "BUILD37", "DRAGEN",
            self.sample_type, name_prep)
        log_base = os.path.join(
            self.scratch_dir, "logs", "{}.{}.{}".format(
                self.pipeline_step_id, name_prep, self.__class__.__name__))
        self.log_file = log_base + ".log"
        self.err = log_base + ".err"
        self.log_dir = os.path.join(self.scratch_dir, "logs")
        self.script = os.path.join(
            self.scratch_dir, "scripts", "{}.{}.{}.sh".format(
                self.pipeline_step_id, name_prep, self.__class__.__name__))
        self.interval_list = os.path.join(
            self.scratch_dir, name_prep + ".interval_list")
        self.scratch_bam = os.path.join(
            self.scratch_dir, name_prep + ".bam")
        self.realn_bam = os.path.join(
            self.scratch_dir, name_prep + ".realn.bam")
        self.recal_table = os.path.join(
            self.scratch_dir, name_prep + ".recal_table")
        self.recal_bam = os.path.join(
            self.scratch_dir, name_prep + ".realn.recal.bam")
        self.recal_bam_index = os.path.join(
            self.scratch_dir, name_prep + ".realn.recal.bai")
        self.gvcf = os.path.join(
            self.scratch_dir, name_prep + ".g.vcf.gz")
        self.gvcf_index = os.path.join(
            self.scratch_dir, name_prep + ".g.vcf.gz.tbi")
        self.vcf = os.path.join(
            self.scratch_dir, name_prep + ".raw.vcf")
        self.snp_vcf = os.path.join(
            self.scratch_dir, name_prep + ".snp.vcf")
        self.indel_vcf = os.path.join(
            self.scratch_dir, name_prep + ".indel.vcf")
        self.snp_recal = os.path.join(
            self.scratch_dir, name_prep + ".snp.recal")
        self.snp_rscript = os.path.join(
            self.scratch_dir, name_prep + ".snp.rscript")
        self.snp_tranches = os.path.join(
            self.scratch_dir, name_prep + ".snp.tranches")
        self.snp_filtered = os.path.join(
            self.scratch_dir, name_prep + ".snp.filtered.vcf")
        self.indel_recal = os.path.join(
            self.scratch_dir, name_prep + ".indel.recal")
        self.indel_rscript = os.path.join(
            self.scratch_dir, name_prep + ".indel.rscript")
        self.indel_tranches = os.path.join(
            self.scratch_dir, name_prep + ".indel.tranches")
        self.indel_filtered = os.path.join(
            self.scratch_dir, name_prep + ".indel.filtered.vcf")
        self.tmp_vcf = os.path.join(
            self.scratch_dir, name_prep + ".tmp.vcf")
        self.final_vcf = os.path.join(
            self.scratch_dir, name_prep + ".analysisReady.vcf")
        self.final_vcf_gz = self.final_vcf + ".gz"
        self.annotated_vcf = os.path.join(
            self.scratch_dir, name_prep + ".analysisReady.annotated.vcf")
        self.annotated_vcf_gz = self.annotated_vcf + ".gz"
        self.annotated_vcf_gz_index = self.annotated_vcf_gz + ".tbi"
        self.phased_vcf = os.path.join(
            self.scratch_dir, name_prep + ".analysisReady.phased.vcf")
        self.phased_vcf_gz = self.phased_vcf + ".gz"
        self.fixed_vcf = os.path.join(
            self.scratch_dir, name_prep + ".analysisReady.fixed.vcf")
        self.original_vcf_gz = os.path.join(
            self.scratch_dir, name_prep +
            ".analysisReady.annotated.original.vcf.gz")
        self.original_vcf_gz_index = self.original_vcf_gz + ".tbi"
        self.alignment_metrics = os.path.join(
            self.scratch_dir, name_prep + ".alignment.metrics.txt")
        self.genome_cov_bed = os.path.join(
            self.scratch_dir, self.sample_name + ".genomecvg.bed")
        self.cov_dir = os.path.join(self.scratch_dir, "cvg_binned")
        self.gq_dir = os.path.join(self.scratch_dir, "gq_binned")
        self.shell_options["record_commands_fn"] = self.script
        self.shell_options["stdout"] = self.log_file
        self.shell_options["stderr"] = self.err
        for d in (os.path.join(self.scratch_dir, "logs"),
                  os.path.join(self.scratch_dir, "scripts")):
            if not os.path.isdir(d):
                os.makedirs(d)
