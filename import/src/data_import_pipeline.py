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
from dragen_globals import *
from db_statements import *
from sys import maxsize

cfg = get_cfg()
CHROMs = OrderedDict([CHROM, x] for x, CHROM in enumerate(
    cfg.get("ref", "CHROMs").split(",")))

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
    def __init__(self, sample_id, pipeline_step_id):
        self.sample_id = sample_id
        self.pipeline_step_id = pipeline_step_id

    def exists(self):
        db = get_connection("dragen")
        try:
            cur = db.cursor()
            cur.execute(GET_STEP_STATUS.format(
                sample_id=self.sample_id,
                pipeline_step_id=self.pipeline_step_id))
            row = cur.fetchone()
            if row:
                if row[0] == 1:
                    return True
                else:
                    return False
            else:
                cur.execute(INSERT_PIPELINE_STEP.format(
                    sample_id=self.sample_id,
                    pipeline_step_id=self.pipeline_step_id))
                db.commit()
                return False
        finally:
            if db.open:
                db.close()

class ParseVCF(SGEJobTask):
    """parse the specified VCF, output text files for the data, and
    import to the database
    """
    vcf = luigi.InputFileParameter(description="the VCF to parse")
    chromosome = luigi.Parameter(description="the chromosome to process")
    sample_id = luigi.NumericalParameter(
        min_value=1, max_value=maxsize, var_type=int,
        description="the sample_id for this sample")
    output_base = luigi.Parameter(
        default=None, description="specify a custom output file name")
    output_directory = luigi.Parameter(
        description="the scratch directory for outputting files")
    sample_name = luigi.Parameter(description="the name of the sample")
    task_name_format = luigi.Parameter(
        significant=False, default="{task_family}_{sample_name}_{chromosome}",
        description="how to format the job name submitted to the cluster")

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
        self.novel_transcripts = self.output_base + ".novel_transcripts.txt"
        self.called_variants = self.output_base + ".calls.txt"
        self.variant_id_vcf = self.output_base + ".variant_id.vcf"
        self.matched_indels = self.output_base + ".matched_indels.txt"
        db = get_connection("dragen")
        try:
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                chromosome=self.chromosome, data_type="Variant Data"))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(CopyDataToScratch)

    def work(self):
        db = get_connection("dragen")
        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                sample_id=self.sample_id,
                pipeline_step_id=self.pipeline_step_id))
            db.commit()
        finally:
            if db.open:
                db.close()
        parse_vcf(
            vcf=self.copied_files_dict["vcf"],
            CHROM=self.chromosome, sample_id=self.sample_id,
            output_base=self.output_base, debug=self.debug)
        for fn in (self.novel_variants, self.novel_transcripts,
                   self.called_variants, self.variant_id_vcf,
                   self.matched_indels):
            if not os.path.isfile(fn):
                raise ValueError("failed running task; {} doesn't exist".format(
                    fn=fn))
        vcf_lines = 0
        vcf_extra_calls = 0
        vcf_tabix = tabix.open(self.copied_files_dict["vcf"])
        for vcf_fields in vcf_tabix.querys(self.chromosome):
            vcf_lines += 1
            vcf_extra_calls += vcf_fields[VCF_COLUMNS_DICT["ALT"]].count(",")
        vcf_variants_count = vcf_lines + vcf_extra_calls
        calls_line_count = get_num_lines_from_vcf(
            self.called_variants, header=False)
        if vcf_variants_count != calls_line_count:
            raise ValueError("incorrect number of variants in calls table")
        variant_id_vcf_line_count = get_num_lines_from_vcf(
            self.variant_id_vcf, header=False)
        if vcf_lines != variant_id_vcf_line_count:
            raise ValueError(
                "incorrect number of lines in variant_id annotated VCF")
        db = get_connection("dragen")
        try:
            cur = db.cursor()
            for table_name, table_file, is_variant_table in (
                ("variant_chr" + self.chromosome, self.novel_variants, True),
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
                    print("error with:")
                    print(load_statement)
                    raise
            cur.execute(GET_NUM_CALLS_FOR_SAMPLE.format(
                CHROM=self.chromosome, sample_id=self.sample_id))
            db_count = cur.fetchone()[0]
            if db_count != vcf_variants_count:
                db.rollback()
                raise ValueError(
                    "incorrect number of calls in the called_variant table")
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                sample_id=self.sample_id,
                pipeline_step_id=self.pipeline_step_id))
            db.commit()
        finally:
            if db.open:
                db.close()

    def output(self):
        return SQLTarget(sample_id=self.sample_id,
                         pipeline_step_id=self.pipeline_step_id)

class LoadGQData(SGEJobTask):
    """Import the GQ binning data for a single chromosome for a sample
    """
    fn = luigi.Parameter(description="the file to import")
    chromosome = luigi.Parameter(description="the chromosome")
    sample_id = luigi.IntParameter(
        description="the sample_id for the sample")
    output_directory = luigi.Parameter(
        description="the scratch directory for outputting files")

    def __init__(self, *args, **kwargs):
        super(LoadGQData, self).__init__(*args, **kwargs)
        db = get_connection("dragen")
        try:
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                chromosome=self.chromosome, data_type="GQ Data"))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
        db = get_connection("dragen")
        cur = db.cursor()
        cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
            sample_id=self.sample_id,
            pipeline_step_id=self.pipeline_step_id))
        db.commit()
        copy(self.fn, self.output_directory)
        temp_fn = os.path.join(self.output_directory, os.path.basename(self.fn))
        try:
            cur.execute(
                INSERT_BIN_STATEMENT.format(
                    data_file=temp_fn, bin_type="GQ",
                    chromosome=self.chromosome, sample_id=self.sample_id))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                sample_id=self.sample_id,
                pipeline_step_id=self.pipeline_step_id))
            db.commit()
        finally:
            if not db.open:
                db.close()

    def output(self):
        return SQLTarget(sample_id=self.sample_id,
                         pipeline_step_id=self.pipeline_step_id)

class LoadDPData(SGEJobTask):
    """Import the DP binning data for a single chromosome for a sample
    """
    fn = luigi.Parameter(description="the file to import")
    chromosome = luigi.Parameter(description="the chromosome")
    sample_id = luigi.IntParameter(
        description="the sample_id for the sample")
    output_directory = luigi.Parameter(
        description="the scratch directory for outputting files")

    def __init__(self, *args, **kwargs):
        super(LoadDPData, self).__init__(*args, **kwargs)
        db = get_connection("dragen")
        try:
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                chromosome=self.chromosome, data_type="DP Data"))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
        db = get_connection("dragen")
        cur = db.cursor()
        cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
            sample_id=self.sample_id,
            pipeline_step_id=self.pipeline_step_id))
        db.commit()
        copy(self.fn, self.output_directory)
        temp_fn = os.path.join(self.output_directory, os.path.basename(self.fn))
        try:
            cur.execute(
                INSERT_BIN_STATEMENT.format(
                    data_file=temp_fn, bin_type="DP",
                    chromosome=self.chromosome, sample_id=self.sample_id))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                sample_id=self.sample_id,
                pipeline_step_id=self.pipeline_step_id))
            db.commit()
        finally:
            if not db.open:
                db.close()

    def output(self):
        return SQLTarget(sample_id=self.sample_id,
                         pipeline_step_id=self.pipeline_step_id)

class ImportSample(luigi.Task):
    """import a sample by requiring() a ParseVCF Task for each chromosome
    """
    vcf = luigi.InputFileParameter(
        default=None,
        description="the VCF to parse (optional - only sample_id is required)")
    sample_id = luigi.NumericalParameter(
        min_value=1, max_value=maxsize, var_type=int,
        description="the sample_id for this sample")
    seqscratch = luigi.ChoiceParameter(
        choices=["09", "10", "11"], default="09",
        description="the seqscratch to use for temporary file creation")
    sample_name = luigi.Parameter(
        default=None,
        description="the name of the sample (doesn't need to be specified")
    output_directory = luigi.Parameter(
        default=None, description="the scratch directory for outputting files "
        "(doesn't need to be specified")
    run_locally = luigi.BoolParameter(
        description="run locally instead of on the cluster")

    def __init__(self, *args, **kwargs):
        super(ImportSample, self).__init__(*args, **kwargs)
        db = get_connection("dragen")
        try:
            cur = db.cursor()
            cur.execute(GET_SAMPLE_INFO.format(sample_id=self.sample_id))
            row = cur.fetchone()
            if row:
                (sample_name, sequencing_type, capture_kit, prep_id) = row
                kwargs["sample_name"] = sample_name
                self.sequencing_type = sequencing_type
                self.prep_id = prep_id
            else:
                raise ValueError("sample_id {sample_id} does not exist".format(
                    sample_id=sample_id))
        finally:
            if db.open:
                db.close()
        kwargs["output_directory"] = cfg.get("pipeline", "scratch_area").format(
            seqscratch=self.seqscratch, sample_name=sample_name,
            sequencing_type=self.sequencing_type.upper())
        self.data_directory = get_data_directory(
            sample_name, self.sequencing_type, capture_kit, self.prep_id)
        if not os.path.isdir(self.data_directory):
            raise ValueError(
                "the data directory {} for the sample does not exist".format(
                self.data_directory))
        db = get_connection("dragen")
        try:
            cur = db.cursor()
            cur.execute(GET_PIPELINE_FINISHED_ID)
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()
        if not self.vcf:
            kwargs["vcf"] = os.path.join(
                self.data_directory, "{sample_name}.{prep_id}.analysisReady."
                "annotated.vcf.gz".format(
                    sample_name=sample_name, prep_id=self.prep_id))
        super(ImportSample, self).__init__(*args, **kwargs)
        if not os.path.isfile(self.vcf):
            raise OSError(
                "could not find the VCF: {vcf}".format(vcf=self.vcf))

    def requires(self):
        return ([self.clone(ParseVCF, chromosome=CHROM, 
            output_base=self.output_directory + CHROM)
            for CHROM in CHROMs.iterkeys()] +
            [self.clone(LoadGQData, fn=os.path.join(
                self.data_directory, "gq_binned",
                "{sample_name}.{prep_id}_gq_binned_10000_chr{chromosome}.txt".format(
                    sample_name=self.sample_name, prep_id=self.prep_id,
                    chromosome=CHROM)), chromosome=CHROM)
             for CHROM in CHROMs.iterkeys()] +
            [self.clone(LoadDPData, fn=os.path.join(
                self.data_directory, "cvg_binned",
                "{sample_name}.{prep_id}_coverage_binned_10000_chr{chromosome}.txt".format(
                    sample_name=self.sample_name, prep_id=self.prep_id,
                    chromosome=CHROM)), chromosome=CHROM)
             for CHROM in CHROMs.iterkeys()])
    
    def run(self):
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
        copy(vcf_out, os.path.join(
            self.data_directory, os.path.basename(vcf_out)))
        copy(vcf_out + ".tbi", os.path.join(
            self.data_directory, os.path.basename(vcf_out) + ".tbi"))
        # copy other stuff/delete scratch directory, etc.
        db = get_connection("dragen")
        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_FINISH_SUBMIT_TIME.format(
                sample_id=self.sample_id,
                pipeline_step_id=self.pipeline_step_id))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                sample_id=self.sample_id,
                pipeline_step_id=self.pipeline_step_id))
            db.commit()
        finally:
            if db.open:
                db.close()
                
    def output(self):
        return SQLTarget(sample_id=self.sample_id,
                         pipeline_step_id=self.pipeline_step_id)

if __name__ == "__main__":
    luigi.run()
