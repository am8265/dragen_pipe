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
import operator
import tabix
import shlex
from shutil import copy
from collections import OrderedDict
from parse_vcf import parse_vcf
from waldb_globals import *
from db_statements import *
import sys
import logging
import Command
from time import sleep

cfg = get_cfg()
# the maximum number of tables to load concurrently
max_tables_to_load = cfg.getint("pipeline", "max_tables_to_load")
table_load_wait_time = cfg.getint("pipeline", "table_load_wait_time")
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
    for data_type in ("Variant", "DP", "GQ"):
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
    cur, max_tables_to_load=max_tables_to_load,
    table_load_wait_time=table_load_wait_time):
    """infinite loop to wait until fewer than the maximum number of tables are
    loading
    """
    while True:
        if get_num_tables_loading(cur) < max_tables_to_load:
            return
        else:
            sleep(table_load_wait_time)

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

    def exists(self):
        for fn in self.targets.itervalues():
            if not os.path.isfile(fn):
                return False
        # verify all files are exactly the same
        with open(os.devnull, "w") as devnull:
            for file_type, fn in self.targets.iteritems():
                p = subprocess.Popen([
                    "diff", "-q", fn, self.originals[file_type]],
                    stdout=devnull, stderr=devnull)
                p.communicate()
                if p.returncode:
                    return False
        return True

    def get_targets(self):
        # just return a dict containing the mapping between file type and new
        # file location
        return self.targets

class CopyDataToScratch(luigi.Task):
    """Copy the VCF and its index to scratch space for processing
    """
    vcf = luigi.InputFileParameter(description="the VCF to copy")
    sample_name = luigi.Parameter(description="the name of the sample")
    output_directory = luigi.Parameter(
        description="the scratch directory for outputting files")

    def __init__(self, *args, **kwargs):
        super(CopyDataToScratch, self).__init__(*args, **kwargs)
        self.vcf_idx = self.vcf + ".tbi"
        self.targets = {}
        self.originals = {}
        for file_type in ("vcf", "vcf_idx"):
            self.targets[file_type] = (
                self.output_directory + os.path.basename(self.__dict__[file_type]))
            self.originals[file_type] = self.__dict__[file_type]

    def run(self):
        for file_type, fn in self.targets.iteritems():
            copy(self.originals[file_type], fn)

    def output(self):
        return CopyDataTarget(targets=self.targets, originals=self.originals)

class SQLTarget(luigi.Target):
    """Target describing verification of the entries in the database
    """
    def __init__(self, prep_id, pipeline_step_id):
        self.prep_id = prep_id
        self.pipeline_step_id = pipeline_step_id

    def exists(self):
        seqdb = get_connection("seqdb")
        try:
            seq_cur = seqdb.cursor()
            seq_cur.execute(GET_STEP_STATUS.format(
                prep_id=self.prep_id,
                pipeline_step_id=self.pipeline_step_id))
            row = seq_cur.fetchone()
            if row:
                if row[0] == 1:
                    return True
                else:
                    return False
            else:
                seq_cur.execute(INSERT_PIPELINE_STEP.format(
                    prep_id=self.prep_id,
                    pipeline_step_id=self.pipeline_step_id))
                seqdb.commit()
                return False
        finally:
            if seqdb.open:
                seqdb.close()

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
        choices=["waldb", "dragen", "waldb4", "waldb1"], default="waldb",
        description="the database to load to")
    min_dp_to_include = luigi.NumericalParameter(
        min_value=0, max_value=sys.maxsize, var_type=int,
        default=cfg.getint("pipeline", "min_dp_to_include"),
        description="ignore variant calls below this read depth")
    dont_load_data = luigi.BoolParameter(
        description="don't actually load any data, used for testing purposes")

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
        self.copied_files_dict = self.input().get_targets()
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
        return self.clone(CopyDataToScratch)

    def work(self):
        if not self.dont_load_data:
            seqdb = get_connection("seqdb")
            try:
                seq_cur = seqdb.cursor()
                seq_cur.execute(GET_TIMES_STEP_RUN.format(
                    prep_id=self.prep_id, pipeline_step_id=self.pipeline_step_id))
                times_run = seq_cur.fetchone()[0] + 1
                seq_cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                    prep_id=self.prep_id, times_run=times_run,
                    pipeline_step_id=self.pipeline_step_id))
                seqdb.commit()
            finally:
                if seqdb.open:
                    seqdb.close()
        level = LOGGING_LEVELS[self.level]
        parse_logger = logging.getLogger("parse_vcf")
        parse_logger.setLevel(level)
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(level)
        formatter = logging.Formatter(cfg.get("logging", "format"))
        handler.setFormatter(formatter)
        parse_logger.addHandler(handler)
        parse_vcf(
            vcf=self.copied_files_dict["vcf"],
            CHROM=self.chromosome, sample_id=self.sample_id,
            database=self.database, min_dp_to_include=self.min_dp_to_include,
            output_base=self.output_base,
            chromosome_length=CHROMs[self.chromosome])
        for fn in (self.novel_variants, self.novel_indels,
                   self.novel_transcripts,
                   self.called_variants, self.variant_id_vcf,
                   self.matched_indels):
            if not os.path.isfile(fn):
                raise ValueError("failed running task; {} doesn't exist".format(
                    fn=fn))
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
                    ("custom_transcript_ids_chr" + self.chromosome,
                     self.novel_transcripts, False),
                    ("called_variant_chr" + self.chromosome, self.called_variants,
                     False),
                    ("matched_indels", self.matched_indels, False)):
                    load_statement = (
                        LOAD_TABLE_REPLACE if is_variant_table else LOAD_TABLE)
                    load_statement = load_statement.format(
                        table_name=table_name, table_file=table_file)
                    try:
                        cur.execute(load_statement)
                    except (MySQLdb.IntegrityError, MySQLdb.OperationalError):
                        logger.error(
                            "error with:\n" + load_statement, exc_info=True)
                        sys.exit(1)
                cur.execute(GET_NUM_CALLS_FOR_SAMPLE.format(
                    CHROM=self.chromosome, sample_id=self.sample_id))
                variant_ids = set()
                for variant_id in cur.fetchall():
                    variant_ids.add(variant_id[0])
                db_count = len(variant_ids)
                if db_count != vcf_variants_count:
                    db.rollback()
                    raise ValueError(
                        "incorrect number of calls in the called_variant table")
                db.commit()
                seq_cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                    prep_id=self.prep_id, pipeline_step_id=self.pipeline_step_id))
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
                prep_id=self.prep_id, pipeline_step_id=self.pipeline_step_id)

class LoadBinData(SGEJobTask):
    """Import the DP/GQ binning data for a single chromosome for a sample
    """
    fn = luigi.Parameter(description="the file to import")
    chromosome = luigi.Parameter(description="the chromosome")
    sample_id = luigi.IntParameter(
        description="the sample_id for the sample")
    output_directory = luigi.Parameter(
        description="the scratch directory for outputting files")
    sample_name = luigi.Parameter(description="the name of the sample")
    prep_id = luigi.IntParameter(
        description="the (pseudo) prep_id for the sample")
    data_type = luigi.ChoiceParameter(
        choices=["DP", "GQ"], description="the type of binned data to load")
    database = luigi.ChoiceParameter(
        choices=["waldb", "dragen", "waldb4", "waldb1"], default="waldb",
        description="the database to load to")
    dont_load_data = luigi.BoolParameter(
        description="don't actually load any data, used for testing purposes")

    def __init__(self, *args, **kwargs):
        super(LoadBinData, self).__init__(*args, **kwargs)
        self.temp_fn = os.path.join(self.output_directory, os.path.basename(self.fn))
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

    def work(self):
        copy(self.fn, self.output_directory)
        if not self.dont_load_data:
            db = get_connection(self.database)
            seqdb = get_connection("seqdb")
            try:
                cur = db.cursor()
                block_until_loading_slot_available(cur)
                seq_cur = seqdb.cursor()
                seq_cur.execute(GET_TIMES_STEP_RUN.format(
                    prep_id=self.prep_id, pipeline_step_id=self.pipeline_step_id))
                times_run = seq_cur.fetchone()[0] + 1
                seq_cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                    times_run=times_run, prep_id=self.prep_id,
                    pipeline_step_id=self.pipeline_step_id))
                seqdb.commit()
                statement = INSERT_BIN_STATEMENT.format(
                    data_file=self.temp_fn, data_type=self.data_type,
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
            return MultipleFilesTarget([[self.temp_fn]])
        else:
            return SQLTarget(
                prep_id=self.prep_id, pipeline_step_id=self.pipeline_step_id)

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
    run_locally = luigi.BoolParameter(
        description="run locally instead of on the cluster")
    level = luigi.ChoiceParameter(
        choices=LOGGING_LEVELS, significant=False,
        default="ERROR", description="the logging level to use")
    dont_remove_tmp_dir_if_failure = luigi.BoolParameter(
        significant=False,
        description="don't remove the tmp dir if there's a failure")
    database = luigi.ChoiceParameter(
        choices=["waldb", "dragen", "waldb4", "waldb1"], default="waldb",
        description="the database to load to")
    min_dp_to_include = luigi.NumericalParameter(
        min_value=0, max_value=sys.maxsize, var_type=int,
        default=cfg.getint("pipeline", "min_dp_to_include"),
        description="ignore variant calls below this read depth")
    dont_load_data = luigi.BoolParameter(
        description="don't actually load any data, used for testing purposes")

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
        kwargs["output_directory"] = cfg.get("pipeline", "scratch_area").format(
            seqscratch=self.seqscratch, sample_name=sample_name,
            prep_id=self.prep_id,
            sequencing_type=self.sequencing_type.upper())
        self.data_directory = get_data_directory(sample_name, self.prep_id)
        try:
            cmd = "ls {} &> /dev/null".format(self.data_directory)
            c = Command.Command(cmd)
            c.run(5)
        except Command.TimeoutException:
            raise ValueError("the data directory {} is inaccessible".format(
                self.data_directory))
        if not os.path.isdir(self.data_directory):
            raise ValueError(
                "the data directory {} for the sample does not exist".format(
                self.data_directory))
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
                self.data_directory, "{sample_name}.{prep_id}.analysisReady."
                "annotated.vcf.gz".format(
                    sample_name=sample_name, prep_id=self.prep_id))
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

    def requires(self):
        return ([self.clone(
            ParseVCF, chromosome=CHROM, 
            output_base=self.output_directory + CHROM)
            for CHROM in CHROMs.iterkeys()] +
            [self.clone(
                LoadBinData,
                fn=os.path.join(
                    self.data_directory, "gq_binned",
                    "{sample_name}.{prep_id}_gq_binned_10000_chr{chromosome}.txt".format(
                        sample_name=self.sample_name, prep_id=self.prep_id,
                        chromosome=CHROM)), chromosome=CHROM, data_type="GQ")
                for CHROM in CHROMs.iterkeys()] +
            [self.clone(
                LoadBinData,
                fn=os.path.join(
                    self.data_directory, "cvg_binned",
                    "{sample_name}.{prep_id}_coverage_binned_1000_chr{chromosome}.txt".format(
                        sample_name=self.sample_name, prep_id=self.prep_id,
                        chromosome=CHROM)), chromosome=CHROM, data_type="DP")
                for CHROM in CHROMs.iterkeys()])
    
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
                seq_cur.execute(INCREMENT_TIMES_STEP_RUN.format(
                    prep_id=self.prep_id, pipeline_step_id=self.pipeline_step_id))
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
        vcf_out = os.path.join(self.output_directory, self.sample_name +
                               ".variant_id.vcf")
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
            cmd = "bgzip " + vcf_out
            p = subprocess.Popen(
                shlex.split(cmd), stdout=devnull, stderr=devnull)
            p.communicate()
            if p.returncode:
                raise subprocess.CalledProcessError(p.returncode, cmd)
            vcf_out += ".gz"
            cmd = "tabix -p vcf " + vcf_out
            p = subprocess.Popen(
                shlex.split(cmd), stdout=devnull, stderr=devnull)
            p.communicate()
            if p.returncode:
                raise subprocess.CalledProcessError(p.returncode, cmd)
        if not self.dont_load_data:
            # only back up/update DB if we're actually loading data
            copy(vcf_out, os.path.join(
                self.data_directory, os.path.basename(vcf_out)))
            copy(vcf_out + ".tbi", os.path.join(
                self.data_directory, os.path.basename(vcf_out) + ".tbi"))
            # copy other stuff/delete scratch directory, etc.
            seqdb = get_connection("seqdb")
            db = get_connection(self.database)
            try:
                seq_cur = seqdb.cursor()
                cur = db.cursor()
                seq_cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                    prep_id=self.prep_id, pipeline_step_id=self.pipeline_step_id))
                cur.execute(SET_SAMPLE_FINISHED.format(
                    sample_id=self.sample_id))
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
                prep_id=self.prep_id, pipeline_step_id=self.pipeline_step_id)

if __name__ == "__main__":
    luigi.run()
