#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Run pipeline to:
    1. Copy the VCF and its index to scratch space
    2. Parse and import data for each chromosome
    3. Archive appropriate data for sample
"""

import luigi
from luigi.contrib.sge import SGEJobTask
import os
import subprocess
import tabix
import shlex
import warnings
import MySQLdb
from shutil import copy
from collections import OrderedDict
from parse_vcf import parse_vcf
from waldb_globals import *
from db_statements import *
import sys
import logging
import Command
import tarfile
from time import sleep
from random import random

warnings.filterwarnings("error", category=MySQLdb.Warning)
cfg = get_cfg()
# the maximum number of tables to load concurrently
MAX_TABLES_TO_LOAD = cfg.getint("pipeline", "max_tables_to_load")
TABLE_LOAD_WAIT_TIME = cfg.getint("pipeline", "table_load_wait_time")
CHROMs = OrderedDict([[chromosome.upper(), int(length)]
                      for chromosome, length in cfg.items("chromosomes")])
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stderr)
formatter = logging.Formatter(cfg.get("logging", "format"))
handler.setFormatter(formatter)
logger.addHandler(handler)

def get_task_list_to_run(cur, prep_id):
    """get the list of all tasks that still need to be run for the sample
    """
    steps_needed = []
    for data_type in ("Variant", "DP"):
        for CHROM in CHROMs:
            step_name = "Chromosome {CHROM} {data_type} Data".format(
                CHROM=CHROM, data_type=data_type)
            cur.execute(STEP_FINISHED.format(
                prep_id=prep_id, step_name=step_name))
            if not cur.fetchone():
                steps_needed.append(step_name)
    return steps_needed

def get_num_tables_loading(cur):
    """return the number of tables that are currently loading in the DB
    """
    num_tables = 0
    cur.execute("SHOW FULL PROCESSLIST")
    for row in cur.fetchall():
        (sql_id, user, host, db, command, time, state, info) = row[:8]
        if db == "WalDB" and info and "LOAD DATA " in info.upper():
            num_tables += 1
    return num_tables

def check_if_sample_importing(cur):
    """check if a sample was presumably stuck previously so we don't try to
    import until the previous one is finished/killed
    """
    return get_num_tables_loading(cur) > 0

def block_until_loading_slot_available(
    cur, max_tables_to_load=MAX_TABLES_TO_LOAD,
    table_load_wait_time=TABLE_LOAD_WAIT_TIME):
    """infinite loop to wait until fewer than the maximum number of tables are
    loading
    """
    while True:
        if get_num_tables_loading(cur) < max_tables_to_load:
            return
        else:
            # sleep a random amount of time, from
            # [(0.5 * table_load_wait_time), (1.5 * table_load_wait_time)]
            base_load_time = table_load_wait_time / 2.0
            wait_time = base_load_time + table_load_wait_time * random()
            sleep(wait_time)

def get_num_lines_from_vcf(vcf, region=None, header=True):
    """return the number of variant lines in the specified VCF, gzipped optional
    """
    count = 0
    if region:
        vcf_tabix = tabix.open(vcf)
        vcf_iter = vcf_tabix.querys(region)
    else:
        vcf_iter = get_fh(vcf)
    if header:
        for line in vcf_iter:
            if line.startswith("#CHROM"):
                break
    for _ in vcf_iter:
        count += 1
    if not region:
        vcf_iter.close()
    return count

class CopyDataTarget(luigi.Target):
    """Target describing the results of copying to scratch space
    """
    def __init__(self, targets, originals):
        self.targets = targets
        self.originals = originals
        assert len(self.targets) == len(self.originals)

    def exists(self):
        for fn in self.targets:
            if not os.path.isfile(fn):
                return False
        # verify all files are exactly the same
        with open(os.devnull, "w") as devnull:
            for target, original in zip(self.targets, self.originals):
                p = subprocess.Popen([
                    "diff", "-q", target, original], stdout=devnull, stderr=devnull)
                p.communicate()
                if p.returncode:
                    return False
        return True

    def get_targets(self):
        # just return a dict containing the mapping between file type and new
        # file location
        return self.targets

class CopyDataToScratch(SGEJobTask):
    """Copy the VCF and its index to scratch space for processing
    """
    originals = luigi.ListParameter(description="the list of files to copy")
    output_directory = luigi.Parameter(
        description="the scratch directory for outputting files")
    sample_name = luigi.Parameter(description="the name of the sample")
    sample_id = luigi.NumericalParameter(
        min_value=1, max_value=sys.maxsize, var_type=int,
        description="the sample_id for this sample")

    def __init__(self, *args, **kwargs):
        kwargs["task_name_format"] = "{task_family}.{sample_name}.{sample_id}"
        super(CopyDataToScratch, self).__init__(*args, **kwargs)
        self.targets = []
        for fn in self.originals:
            self.targets.append(os.path.join(
                self.output_directory, os.path.basename(fn)))

    def work(self):
        for original, target in zip(self.originals, self.targets):
            copy(original, target)

    def output(self):
        return CopyDataTarget(targets=self.targets, originals=self.originals)

class ExtractTarball(SGEJobTask):
    """Extract a tarball to the output directory
    """
    tarball = luigi.Parameter(
        description="the path to the tarball to copy and extract")
    fn_format = luigi.Parameter(
        description="the format of the file names in the tarball")
    output_directory = luigi.Parameter(
        description="the scratch directory for outputting files")
    sample_name = luigi.Parameter(description="the name of the sample")
    prep_id = luigi.IntParameter(
        description="the (pseudo) prep_id for the sample")
    sample_id = luigi.NumericalParameter(
        min_value=1, max_value=sys.maxsize, var_type=int,
        description="the sample_id for this sample")

    def __init__(self, *args, **kwargs):
        kwargs["task_name_format"] = "{task_family}.{sample_name}.{sample_id}"
        super(ExtractTarball, self).__init__(*args, **kwargs)

    def requires(self):
        return self.clone( CopyDataToScratch, originals=[self.tarball] )

    def work(self):
        new_tarball = os.path.join(
            self.output_directory, os.path.basename(self.tarball))
        with tarfile.open(new_tarball, "r:gz") as tar:
            tar.extractall(self.output_directory)

    def output(self):
        return MultipleFilesTarget(
            [[os.path.join(self.output_directory, self.fn_format.format(
                chromosome=chromosome, **self.__dict__)) for chromosome in
                CHROMs]])

class ParseVCF(SGEJobTask):
    """parse the specified VCF, output text files for the data, and
    import to the database
    """
    vcf = luigi.InputFileParameter(description="the VCF to parse")
    chromosome = luigi.Parameter(description="the chromosome to process")
    sample_id = luigi.NumericalParameter(
        min_value=1, max_value=sys.maxsize, var_type=int,
        description="the sample_id for this sample")
    output_base = luigi.Parameter(
        default=None, description="specify a custom output file name")
    output_directory = luigi.Parameter(
        description="the scratch directory for outputting files")
    sample_name = luigi.Parameter(description="the name of the sample")
    prep_id = luigi.IntParameter(
        description="the (pseudo) prep_id for the sample")
    level = luigi.ChoiceParameter(
        choices=LOGGING_LEVELS, significant=False,
        default="DEBUG", description="the logging level to use")
    database = luigi.ChoiceParameter(
        choices=["waldb4", "waldb6"],
        default="waldb6", description="the database to load to")
    min_dp_to_include = luigi.NumericalParameter(
        min_value=0, max_value=sys.maxsize, var_type=int,
        default=cfg.getint("pipeline", "min_dp_to_include"),
        description="ignore variant calls below this read depth")
    dont_load_data = luigi.BoolParameter(
        description="don't actually load any data, used for testing purposes")
    priority = 1 # tell the scheduler to run these first
    version = luigi.Parameter(
        description="the version of the pipeline run")

    def __init__(self, *args, **kwargs):
        super(ParseVCF, self).__init__(*args, **kwargs)
        if self.chromosome not in CHROMs:
            raise ValueError("invalid chromosome specified: {chromosome}".format(
                chromosome=self.chromosome))
        if self.output_base:
            if not os.path.isdir(os.path.dirname(self.output_base)):
                os.makedirs(os.path.dirname(self.output_base))
        else:
            kwargs["output_base"] = os.path.splitext(self.vcf)[0]
            super(ParseVCF, self).__init__(*args, **kwargs)
        self.copied_files  = self.input().get_targets()
        self.novel_variants = self.output_base + ".novel_variants.txt"
        self.novel_indels = self.output_base + ".novel_indels.txt"
        self.novel_transcripts = self.output_base + ".novel_transcripts.txt"
        self.called_variants = self.output_base + ".calls.txt"
        self.variant_id_vcf = self.output_base + ".variant_id.vcf"
        self.matched_indels = self.output_base + ".matched_indels.txt"
        if not self.dont_load_data:
            seqdb = get_connection("seqdb")
            try:
                seq_cur = seqdb.cursor()
                seq_cur.execute(GET_DATA_LOADED_PIPELINE_STEP_ID.format(
                    chromosome=self.chromosome, data_type="Variant Data"))
                self.pipeline_step_id = seq_cur.fetchone()[0]
            finally:
                if seqdb.open:
                    seqdb.close()

    def requires(self):
        return self.clone( CopyDataToScratch, originals=[self.vcf, self.vcf + ".tbi"] )

    def work(self):
        if not self.dont_load_data:
            seqdb = get_connection("seqdb")
            try:
                seq_cur = seqdb.cursor()
                seq_cur.execute(GET_TIMES_STEP_RUN.format(
                    prep_id=self.prep_id, pipeline_step_id=self.pipeline_step_id))
                row = seq_cur.fetchone()
                times_ran = row[0] if row else 0
                seq_cur.execute(BEGIN_STEP.format(
                    pseudo_prepid=self.prep_id,
                    pipeline_step_id=self.pipeline_step_id,
                    version=self.version, times_ran=times_ran + 1))
                seqdb.commit()
            finally:
                if seqdb.open:
                    seqdb.close()
        level = LOGGING_LEVELS[self.level]
        parse_logger = logging.getLogger("parse_vcf")
        parse_logger.setLevel(level)
        parse_vcf_handler = logging.StreamHandler(sys.stdout)
        parse_vcf_handler.setLevel(level)
        handler.setFormatter(formatter)
        parse_logger.addHandler(handler)
        parse_vcf(
            vcf=self.copied_files[0],
            CHROM=self.chromosome, sample_id=self.sample_id,
            database=self.database, min_dp_to_include=self.min_dp_to_include,
            output_base=self.output_base)
        for fn in (self.novel_variants, self.novel_indels,
                   self.novel_transcripts,
                   self.called_variants, self.variant_id_vcf,
                   self.matched_indels):
            if not os.path.isfile(fn):
                import ipdb
                ipdb.set_trace()
                ipdb.pm()
                raise ValueError("failed running task; {} doesn't exist".format(fn=fn))
        variants = set()
        with open(self.variant_id_vcf) as vcf:
            for line in vcf:
                fields = VCF_fields_dict(line.strip("\n").split("\t"))
                dp = int(dict(zip(fields["FORMAT"].split(":"),
                                  fields["call"].split(":")))["DP"])
                if dp < self.min_dp_to_include:
                    # make sure not to count variants with depth less than 3
                    continue
                info = fields["INFO"]
                if info.startswith("VariantID="):
                    for variant_id in info.split(";")[0].split("=")[1].split(","):
                        variants.add(int(variant_id))
                else:
                    raise ValueError("error with format of variant_id VCF at "
                                     "{CHROM}-{POS}-{REF}-{ALT}".format(**fields))
        vcf_variants_count = len(variants)
        calls_line_count = get_num_lines_from_vcf(
            self.called_variants, header=False)
        if vcf_variants_count != calls_line_count:
            raise ValueError("incorrect number of variants in calls table: "
                             "{vcf_count} in VCF, {table_count} in table".format(
                                 vcf_count=vcf_variants_count,
                                 table_count=calls_line_count))
        if not self.dont_load_data:
            db = get_connection(self.database)
            seqdb = get_connection("seqdb")
            try:
                cur = db.cursor()
                seq_cur = seqdb.cursor()
                block_until_loading_slot_available(cur)
                for table_name, table_file, is_variant_table in (
                    ("variant_chr" + self.chromosome, self.novel_variants, True),
                    ("indel_chr" + self.chromosome, self.novel_indels, False),
                    ("custom_transcript_ids_chr" + self.chromosome, self.novel_transcripts, False),
                    ("called_variant_chr" + self.chromosome, self.called_variants, False),
                    ("matched_indels", self.matched_indels, False)):
                    load_statement = ( LOAD_TABLE_REPLACE if is_variant_table else LOAD_TABLE)
                    load_statement = load_statement.format( table_name=table_name, table_file=table_file )
                    try:
                        cur.execute(load_statement)
                    except (MySQLdb.IntegrityError, MySQLdb.OperationalError):
                        logger.error(
                            "error with:\n" + load_statement, exc_info=True)
                        sys.exit(1)
                cur.execute(GET_NUM_CALLS_FOR_SAMPLE.format( CHROM=self.chromosome, sample_id=self.sample_id) )
                variant_ids = set()
                for variant_id in cur.fetchall():
                    variant_ids.add(variant_id[0])
                db_count = len(variant_ids)

                if db_count != vcf_variants_count:
                    db.rollback()
                    raise ValueError( "incorrect number of calls in the called_variant table" )
                db.commit()
                seq_cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format( prep_id=self.prep_id, pipeline_step_id=self.pipeline_step_id) )
                seqdb.commit()

            finally:
                if db.open:
                    db.close()
                if seqdb.open:
                    seqdb.close()

    def output(self):
        if self.dont_load_data:
            return MultipleFilesTarget(
                [[self.novel_variants, self.novel_indels, self.novel_transcripts,
                  self.called_variants, self.variant_id_vcf, self.matched_indels]])
        else:
            return SQLTarget(
                pseudo_prepid=self.prep_id, pipeline_step_id=self.pipeline_step_id)

class LoadBinData(SGEJobTask):
    """Import the DP/GQ binning data for a single chromosome for a sample
    """
    fn = luigi.Parameter(default=None, description="the file to import")
    chromosome = luigi.Parameter(description="the chromosome")
    sample_id = luigi.IntParameter(
        description="the sample_id for the sample")
    output_directory = luigi.Parameter(
        description="the scratch directory for outputting files")
    data_directory = luigi.Parameter(
        default=None, description="the sample's data directory "
        "(doesn't need to be specified)")
    sample_name = luigi.Parameter(description="the name of the sample")
    prep_id = luigi.IntParameter(
        description="the (pseudo) prep_id for the sample")
    data_type = luigi.ChoiceParameter(
        choices=["DP", "GQ"], description="the type of binned data to load")
    database = luigi.ChoiceParameter(
        choices=["waldb_master","waldb", "dragen", "waldb4", "waldb1", "waldb6"],
        default="waldb6", description="the database to load to")
    dont_load_data = luigi.BoolParameter(
        description="don't actually load any data, used for testing purposes")
    version = luigi.Parameter(
        description="the version of the pipeline run")

    def __init__(self, *args, **kwargs):
        super(LoadBinData, self).__init__(*args, **kwargs)
        if self.data_type == "DP":
            self.fn_format = (
                "{sample_name}.{prep_id}_coverage_binned_1000_chr{chromosome}.txt")
        else:
            raise NotImplementedError(
                "This type of bin {} is not supported!".format(self.data_type))
        self.fn = os.path.join(
            self.output_directory, self.fn_format.format(
                sample_name=self.sample_name, prep_id=self.prep_id,
                chromosome=self.chromosome))
        if not self.dont_load_data:
            seqdb = get_connection("seqdb")
            try:
                seq_cur = seqdb.cursor()
                seq_cur.execute(GET_DATA_LOADED_PIPELINE_STEP_ID.format(
                    chromosome=self.chromosome, data_type="{data_type} Data".format(
                    data_type=self.data_type)))
                self.pipeline_step_id = seq_cur.fetchone()[0]
            finally:
                if seqdb.open:
                    seqdb.close()

    def requires(self):
        return self.clone( ExtractTarball, fn_format=self.fn_format, tarball=os.path.join(self.data_directory, "coverage.tar.gz") )

    def work(self):
        if not self.dont_load_data:
            db = get_connection(self.database)
            seqdb = get_connection("seqdb")
            try:
                cur = db.cursor()
                block_until_loading_slot_available(cur)
                seq_cur = seqdb.cursor()
                seq_cur.execute(GET_TIMES_STEP_RUN.format(
                    prep_id=self.prep_id, pipeline_step_id=self.pipeline_step_id))
                row = seq_cur.fetchone()
                times_ran = row[0] if row else 0
                seq_cur.execute(BEGIN_STEP.format(
                    pseudo_prepid=self.prep_id,
                    pipeline_step_id=self.pipeline_step_id,
                    version=self.version, times_ran=times_ran + 1))
                seqdb.commit()
                statement = INSERT_BIN_STATEMENT.format(
                    data_file=self.fn, data_type=self.data_type,
                    chromosome=self.chromosome, sample_id=self.sample_id)
                cur.execute(statement)
                db.commit()
                seq_cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                    prep_id=self.prep_id, pipeline_step_id=self.pipeline_step_id))
                seqdb.commit()
            except MySQLdb.InternalError:
                logger.error("{statement} failed".format(statement=statement))
            finally:
                if db.open:
                    db.close()
                if seqdb.open:
                    seqdb.close()

    def output(self):
        if self.dont_load_data:
            return MultipleFilesTarget([[self.fn]])
        else:
            return SQLTarget(
                pseudo_prepid=self.prep_id, pipeline_step_id=self.pipeline_step_id)

class ImportSample(luigi.Task):
    """import a sample by requiring() a ParseVCF Task for each chromosome
    """
    vcf = luigi.InputFileParameter(
        default=None,
        description="the VCF to parse (optional - only sample_id is required)")
    sample_id = luigi.NumericalParameter(
        min_value=1, max_value=sys.maxsize, var_type=int,
        description="the sample_id for this sample")
    prep_id = luigi.IntParameter(
        default=None,
        description="the (pseudo) prep_id for the sample")
    seqscratch = luigi.ChoiceParameter(
        choices=["09", "10", "11", "_ssd"], default="_ssd",
        description="the seqscratch to use for temporary file creation")
    sample_name = luigi.Parameter(
        default=None,
        description="the name of the sample (doesn't need to be specified")
    output_directory = luigi.Parameter(
        default=None, description="the scratch directory for outputting files "
        "(doesn't need to be specified")
    data_directory = luigi.Parameter(
        default=None, description="the sample's data directory "
        "(doesn't need to be specified)")
    run_locally = luigi.BoolParameter(
        description="run locally instead of on the cluster")
    level = luigi.ChoiceParameter(
        choices=LOGGING_LEVELS, significant=False,
        default="ERROR", description="the logging level to use")
    dont_remove_tmp_dir_if_failure = luigi.BoolParameter(
        significant=False,
        description="don't remove the tmp dir if there's a failure")
    database = luigi.ChoiceParameter(
        choices=["waldb_master","waldb", "dragen", "waldb4", "waldb1", "waldb6"],
        default="waldb6", description="the database to load to")
    min_dp_to_include = luigi.NumericalParameter(
        min_value=0, max_value=sys.maxsize, var_type=int,
        default=cfg.getint("pipeline", "min_dp_to_include"),
        description="ignore variant calls below this read depth")
    dont_load_data = luigi.BoolParameter(
        description="don't actually load any data, used for testing purposes")
    override_directory_check = luigi.BoolParameter(
        description="don't check for existence of sample directory "
        "(will assume it exists and copy files; can be helpful if ls hangs)")
    version = luigi.Parameter(
        default=None, description="the version of the pipeline run "
        "(doesn't need to be specified)")

    def __init__(self, *args, **kwargs):
        super(ImportSample, self).__init__(*args, **kwargs)
        db = get_connection(self.database)
        try:
            cur = db.cursor()

            cur.execute(GET_SAMPLE_INFO.format(sample_id=self.sample_id))
            row = cur.fetchone()
            if row:
                (sample_name, sequencing_type, capture_kit, prep_id) = row
                kwargs["sample_name"] = sample_name
                self.sample_name = sample_name
                self.sequencing_type = sequencing_type
                self.prep_id = prep_id
                if "prep_id" not in kwargs:
                    kwargs["prep_id"] = prep_id
            else:
                raise ValueError("sample_id {sample_id} does not exist".format(
                    sample_id=self.sample_id))
        finally:
            if db.open:
                db.close()

        try:
            seqdb = get_connection("seqdb")
            seq_cur = seqdb.cursor()
            seq_cur.execute('select is_merged from dragen_sample_metadata where pseudo_prepid = {}'.format(self.prep_id))
            ##### horrid but in a hurry and just don't care!?!
            if seq_cur.fetchone()[0]!=41:
                raise ValueError("please don't run me directly ({})!?! ".format(self.prep_id))
            else:
                # should lock the sample itself in waldb to a user/host/pid?!?
                # print('YAY')
                ##### horrid but in a hurry and just don't care!?!
                seq_cur.execute('update dragen_sample_metadata set is_merged = 42 where pseudo_prepid = {} '.format(self.prep_id))
                seqdb.commit()

        finally:
            if db.open:
                db.close()

        kwargs["output_directory"] = cfg.get("pipeline", "scratch_area").format(
            seqscratch=self.seqscratch, sample_name=sample_name,
            prep_id=self.prep_id,
            sequencing_type=self.sequencing_type.upper())
        data_directory = get_data_directory(sample_name, self.prep_id)
        kwargs["data_directory"] = data_directory
        if not self.override_directory_check:
            # this can optionally be disabled if ls hangs but copying files is
            # still working
            try:
                cmd = "ls -d {} &> /dev/null".format(data_directory)
                c = Command.Command(cmd, shell=True)
                c.start()
                #### way too short for something as brute force as this on an object store...
                c.join(15)
                # c.join(5)
            except Command.TimeoutException:
                raise ValueError("the data directory {} is inaccessible".format(
                    data_directory))
            if not os.path.isdir(data_directory):
                raise ValueError(
                    "the data directory {} for the sample does not exist".format(
                    data_directory))
        db = get_connection(self.database)
        seqdb = get_connection("seqdb")
        try:
            cur = db.cursor()
            seq_cur = seqdb.cursor()
            while check_if_sample_importing(cur):
                logger.warning(
                    "Data still loading for other task...sleeping 10 seconds")
                sleep(10)
            logger.debug("Loading {sample_name}:{sequencing_type}:{capture_kit}:"
                         "{prep_id}".format(capture_kit=capture_kit, **self.__dict__))
            seqdb = get_connection("seqdb")
            seq_cur = seqdb.cursor()
            if not self.dont_load_data:
                logger.debug("The following still need to be loaded: {}".format(
                    get_task_list_to_run(seq_cur, self.prep_id)))
            seq_cur.execute(GET_PIPELINE_FINISHED_ID)
            self.pipeline_step_id = seq_cur.fetchone()[0]
        finally:
            if db.open:
                db.close()
            if seqdb.open:
                seqdb.close()
        if not self.vcf:
            kwargs["vcf"] = os.path.join(
                data_directory, "{sample_name}.{prep_id}.analysisReady."
                "annotated.vcf.gz".format(
                    sample_name=sample_name, prep_id=self.prep_id))
        if not self.version:
            kwargs["version"] = get_pipeline_version()
        super(ImportSample, self).__init__(*args, **kwargs)
        if os.path.isfile(self.vcf):
            if self.sequencing_type == "exome" and os.path.getsize(self.vcf) > 248000000:
                db = get_connection(self.database)
                try:
                    cur = db.cursor()
                    cur.execute(
                        "UPDATE sample SET sample_failure = 1 WHERE sample_id = "
                        "{sample_id}".format(sample_id=self.sample_id))
                    db.commit()
                    db.close()
                    raise ValueError(
                        "{sample_name} appears to be a genome sample".
                        format(sample_name=self.sample_name))
                finally:
                    if db.open:
                        db.close()
        else:
            raise OSError(
                "could not find the VCF: {vcf}".format(vcf=self.vcf))

    def requires(self): # list of parse & loadbin 
        return ( 
                    [ self.clone( ParseVCF, chromosome=CHROM,  output_base=self.output_directory + CHROM) for CHROM in CHROMs.iterkeys() ] 
                +
                    [ self.clone( LoadBinData, chromosome=CHROM, data_type="DP") for CHROM in CHROMs.iterkeys() ]
        )
    
    def run(self):
        # verify that each VariantID annotated chromosome's VCF is present,
        # otherwise re-create
        variant_id_vcfs_missing = []
        for CHROM in CHROMs.iterkeys():
            if not os.path.isfile(os.path.join(
                self.output_directory, CHROM + ".variant_id.vcf")):
                variant_id_vcfs_missing.append(CHROM)
        if variant_id_vcfs_missing:
            yield [self.clone(ParseVCF, chromosome=CHROM,
                              output_base=self.output_directory + CHROM,
                              dont_load_data=True) for CHROM in
                   variant_id_vcfs_missing]
        if not self.dont_load_data:
            seqdb = get_connection("seqdb")
            try:
                seq_cur = seqdb.cursor()
                seq_cur.execute(GET_TIMES_STEP_RUN.format(
                    prep_id=self.prep_id, pipeline_step_id=self.pipeline_step_id))
                row = seq_cur.fetchone()
                times_ran = row[0] if row else 0
                seq_cur.execute(BEGIN_STEP.format(
                    pseudo_prepid=self.prep_id,
                    pipeline_step_id=self.pipeline_step_id,
                    version=self.version, times_ran=times_ran + 1))
                seqdb.commit()
            finally:
                if seqdb.open:
                    seqdb.close()
        variant_id_header = cfg.get("pipeline", "variant_id_header") + "\n"
        header = []
        info_encountered = False
        info_output = False
        with get_fh(self.vcf) as vcf_fh:
            for line in vcf_fh:
                if line.startswith("#CHROM"):
                    break
                if line.startswith("##INFO=<"):
                    info_encountered = True
                    header.append(line)
                else:
                    if info_encountered and not info_output:
                        header.append(variant_id_header)
                        info_output = True
                    header.append(line)
        if not info_output:
            header.append(variant_id_header)
        header.append(line)
        vcf_out = os.path.join(
            self.output_directory,
            "{}.{}.variant_id.vcf".format(self.sample_name, self.prep_id))
        with open(vcf_out, "w") as vcf_out_fh:
            for line in header:
                vcf_out_fh.write(line)
            for CHROM in CHROMs.iterkeys():
                with open(os.path.join(
                    self.output_directory, CHROM + ".variant_id.vcf")) as vcf_fh:
                    for line in vcf_fh:
                        vcf_out_fh.write(line)
        with open(os.devnull, "w") as devnull:
            for fn in (vcf_out + ".gz", vcf_out + ".gz.tbi"):
                # attempting to bgzip if the file already exists will cause this
                # to hang indefinitely
                if os.path.isfile(fn):
                    os.remove(fn)
            ##### bored, bored, bored!?!
            cmd = "/nfs/goldstein/software/bin/bgzip " + vcf_out
            print("running : '{}'".format(cmd))
            p = subprocess.Popen( shlex.split(cmd), stdout=devnull, stderr=devnull )
            p.communicate()
            if p.returncode:
                raise subprocess.CalledProcessError(p.returncode, cmd)
            vcf_out += ".gz"
            ##### bored, bored, bored!?!
            cmd = "/nfs/goldstein/software/bin/tabix -p vcf " + vcf_out
            print("running : '{}'".format(cmd))
            p = subprocess.Popen( shlex.split(cmd), stdout=devnull, stderr=devnull )
            p.communicate()
            if p.returncode:
                raise subprocess.CalledProcessError(p.returncode, cmd)
        if not self.dont_load_data:
            ##### currently this is have to run 'after' archive protection which means this is a bit of a PIA and is simpler to move elsewhere or not bother?!?
            # only back up/update DB if we're actually loading data
            # copy(vcf_out, os.path.join( self.data_directory, os.path.basename(vcf_out)))
            # copy(vcf_out + ".tbi", os.path.join( self.data_directory, os.path.basename(vcf_out) + ".tbi"))
            # copy other stuff/delete scratch directory, etc.
            seqdb = get_connection("seqdb")
            db = get_connection(self.database)
            try:
                seq_cur = seqdb.cursor()
                cur = db.cursor()
                seq_cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format( prep_id=self.prep_id, pipeline_step_id=self.pipeline_step_id) )
                seq_cur.execute(UPDATE_PREP_STATUS.format( pseudo_prepid=self.prep_id, status="In DragenDB") )
                ##### horrid but in a hurry and just don't care!?!
                seq_cur.execute('update dragen_sample_metadata set is_merged = 43 where pseudo_prepid = {} '.format(self.prep_id))
                cur.execute(SET_SAMPLE_FINISHED.format( sample_id=self.sample_id) )
                seqdb.commit()
                db.commit()
            finally:
                if seqdb.open:
                    seqdb.close()
                if db.open:
                    db.close()
                
    def output(self):
        if self.dont_load_data:
            # just check for existence of all files
            return MultipleFilesTarget(
                [target.get_targets() for target in self.input()])
        else:
            # check final DB status
            return SQLTarget(
                pseudo_prepid=self.prep_id, pipeline_step_id=self.pipeline_step_id)

if __name__ == "__main__":
    sys.exit(luigi.run())
