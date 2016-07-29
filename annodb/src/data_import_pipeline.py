#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Run pipeline to:
    1. Copy VCF, gVCF, and pileup to scratch space
    2. Parse and import data for each chromosome
    3. Archive appropriate data for sample
"""

import luigi
from luigi.contrib.sge import SGEJobTask
import os
import subprocess
import operator
import tabix
from shutil import copy
from collections import OrderedDict
from parse_vcf import parse_vcf
from dragen_globals import *
from db_statements import *

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
    """Copy the VCF, gVCF, pileup, and their indexes to scratch space for
    processing
    """
    vcf = luigi.InputFileParameter(description="the VCF to copy")
    gvcf = luigi.InputFileParameter(description="the gVCF to copy")
    pileup = luigi.InputFileParameter(description="the pileup to copy")
    seqscratch = luigi.ChoiceParameter(
        choices=["09", "10", "11"], default="09",
        description="the seqscratch to use for temporary file creation")
    sample_name = luigi.Parameter(description="the name of the sample")
    sequencing_type = luigi.ChoiceParameter(
        choices=["exome", "genome", "custom_capture", "merged"],
        description="the type of sequencing performed on this sample")

    def __init__(self, *args, **kwargs):
        super(CopyDataToScratch, self).__init__(*args, **kwargs)
        self.vcf_idx = self.vcf + ".tbi"
        self.gvcf_idx = self.gvcf + ".tbi"
        self.pileup_idx = self.pileup + ".tbi"
        self.output_directory = cfg.get("pipeline", "scratch_area").format(
            seqscratch=self.seqscratch, sample_name=self.sample_name,
            sequencing_type=self.sequencing_type.upper())
        self.targets = {}
        self.originals = {}
        for file_type in (
            "vcf", "gvcf", "pileup", "vcf_idx", "gvcf_idx", "pileup_idx"):
            self.targets[file_type] = (
                self.output_directory + os.path.basename(self.__dict__[file_type]))
            self.originals[file_type] = self.__dict__[file_type]

    def run(self):
        for file_type, fn in self.targets.iteritems():
            copy(self.originals[file_type], fn)

    def output(self):
        return CopyDataTarget(targets=self.targets, originals=self.originals)

class ParsedVCFTarget(luigi.Target):
    """Target describing verification of the entries in the database
    """
    def __init__(self, chromosome, sample_id):
        self.chromosome = chromosome
        self.sample_id = sample_id

    def exists(self):
        db = get_connection("dragen")
        try:
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                chromosome=self.chromosome))
            pipeline_step_id = cur.fetchone()[0]
            cur.execute(GET_STEP_STATUS.format(
                sample_id=self.sample_id,
                pipeline_step_id=pipeline_step_id))
            row = cur.fetchone()
            if row:
                if row[0] == 1:
                    return True
                else:
                    return False
            else:
                cur.execute(INSERT_PIPELINE_STEP.format(
                    sample_id=self.sample_id,
                    pipeline_step_id=pipeline_step_id))
                db.commit()
                return False
        finally:
            if db.open:
                db.close()

class ParseVCF(SGEJobTask):
    """parse the specified VCF and gVCF, output text files for the data, and
    import to the database
    """
    vcf = luigi.InputFileParameter(description="the VCF to parse")
    gvcf = luigi.InputFileParameter(description="the gVCF to parse")
    pileup = luigi.InputFileParameter(description="the pileup to parse")
    chromosome = luigi.Parameter(description="the chromosome to process")
    sample_id = luigi.NumericalParameter(
        left_op=operator.le, description="the sample_id for this sample")
    output_base = luigi.Parameter(
        default=None, description="specify a custom output file name")
    seqscratch = luigi.ChoiceParameter(
        choices=["09", "10", "11"], default="09",
        description="the seqscratch to use for temporary file creation")
    sample_name = luigi.Parameter(description="the name of the sample")
    sequencing_type = luigi.ChoiceParameter(
        choices=["exome", "genome", "custom_capture", "merged"],
        description="the type of sequencing performed on this sample")
    task_name_format = luigi.Parameter(
        significant=False, default="{task_family}_{sample_name}_{chromosome}",
        description="how to format the job name submitted to the cluster")
    debug = luigi.BoolParameter(
        significant=False,
        description="output debug timing messages to job.err")

    def __init__(self, *args, **kwargs):
        super(ParseVCF, self).__init__(*args, **kwargs)
        if self.chromosome not in CHROMs:
            raise ValueError("invalid chromosome specified: {chromosome}".format(
                chromosome=self.chromosome))
        if self.output_base:
            if not os.path.isdir(os.path.dirname(self.output_base)):
                os.makedirs(os.path.dirname(self.output_base))
        else:
            self.output_base = os.path.splitext(self.vcf)[0]
        self.copied_files_dict = self.input().get_targets()
        self.novel_variants = self.output_base + ".novel_variants.txt"
        self.called_variants = self.output_base + ".calls.txt"
        self.variant_id_vcf = self.output_base + ".variant_id.vcf"
        self.matched_indels = self.output_base + ".matched_indels.txt"

    def requires(self):
        return self.clone(CopyDataToScratch)

    def work(self):
        db = get_connection("dragen")
        try:
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(chromosome=self.chromosome))
            pipeline_step_id = cur.fetchone()[0]
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                sample_id=self.sample_id, pipeline_step_id=pipeline_step_id))
            db.commit()
        finally:
            if db.open:
                db.close()
        parse_vcf(
            vcf=self.copied_files_dict["vcf"],
            gvcf=self.copied_files_dict["gvcf"],
            pileup=self.copied_files_dict["pileup"],
            CHROM=self.chromosome, sample_id=self.sample_id,
            output_base=self.output_base, debug=self.debug)
        for fn in (self.novel_variants, self.called_variants,
                   self.variant_id_vcf, self.matched_indels):
            if not os.path.isfile(fn):
                raise ValueError("failed running task; {} doesn't exist".format(
                    fn=fn))
        vcf_lines = 0
        vcf_extra_calls = 0
        vcf_tabix = tabix.open(self.vcf)
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
            for table_name, table_file in (
                ("variant_chr" + self.chromosome, self.novel_variants),
                ("called_variant_chr" + self.chromosome, self.called_variants),
                ("matched_indels", self.matched_indels)):
                cur.execute(LOAD_TABLE.format(
                    table_name=table_name, table_file=table_file))
            cur.execute(GET_NUM_CALLS_FOR_SAMPLE.format(
                CHROM=self.chromosome, sample_id=self.sample_id))
            db_count = cur.fetchone()[0]
            if db_count != vcf_variants_count:
                db.rollback()
                raise ValueError(
                    "incorrect number of calls in the called_variant table")
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                sample_id=self.sample_id, pipeline_step_id=pipeline_step_id))
            db.commit()
        finally:
            if db.open:
                db.close()

    def output(self):
        return ParsedVCFTarget(
            chromosome=self.chromosome, sample_id=self.sample_id)

class ImportSample(luigi.WrapperTask):
    """import a sample by requiring() a ParseVCF Task for each chromosome
    """
    vcf = luigi.InputFileParameter(
        default=None,
        description="the VCF to parse (optional - only sample_id is required)")
    gvcf = luigi.InputFileParameter(
        default=None,
        description="the gVCF to parse (optional - only sample_id is required)")
    pileup = luigi.InputFileParameter(
        default=None,
        description="the pileup to parse (optional - only sample_id is required)")
    sample_id = luigi.NumericalParameter(
        left_op=operator.le, description="the sample_id for this sample")
    seqscratch = luigi.ChoiceParameter(
        choices=["09", "10", "11"], default="09",
        description="the seqscratch to use for temporary file creation")

    def __init__(self, *args, **kwargs):
        super(ImportSample, self).__init__(*args, **kwargs)
        db = get_connection("dragen")
        try:
            cur = db.cursor()
            cur.execute(GET_SAMPLE_INFO.format(sample_id=self.sample_id))
            row = cur.fetchone()
            if row:
                (self.sample_name, self.sequencing_type,
                 self.capture_kit, self.prep_id) = row
            else:
                raise ValueError("sample_id {sample_id} does not exist".format(
                    sample_id=self.sample_id))
        finally:
            if db.open:
                db.close()
        self.output_directory = cfg.get("pipeline", "scratch_area").format(
            seqscratch=self.seqscratch, sample_name=self.sample_name,
            sequencing_type=self.sequencing_type.upper())
        for fn in (self.vcf, self.gvcf, self.pileup):
            if not fn:
                db = get_connection("seqdb")
                try:
                    cur = db.cursor()
                    cur.execute(GET_DATA_DIRECTORY_FOR_SAMPLE.format(
                        sample_name=self.sample_name,
                        sample_type=self.sequencing_type,
                        capture_kit=self.capture_kit, prep_id=self.prep_id))
                    row = cur.fetchone()
                    if row:
                        data_dir = row[0]
                        self.vcf = os.path.join(
                            data_dir, "{sample_name}.hard-filtered.ann.vcf.gz".
                            format(sample_name=self.sample_name))
                        self.gvcf = os.path.join(
                            data_dir, "{sample_name}.hard-filtered.gvcf.gz".
                            format(sample_name=self.sample_name))
                        self.pileup = os.path.join(
                            data_dir, "{sample_name}.genomecvg.bed.gz".format(
                                sample_name=self.sample_name))
                        for file_type in ("vcf", "gvcf", "pileup"):
                            if not os.path.isfile(self.__dict__[file_type]):
                                raise OSError(
                                    "could not find {file_type}: {fn}".format(
                                        file_type=file_type,
                                        fn=self.__dict__[file_type]))
                    else:
                        raise ValueError(
                            "could not find data directory for sample")
                finally:
                    if db.open:
                        db.close()

    def requires(self):
        return [ParseVCF(
            vcf=self.vcf, gvcf=self.gvcf, pileup=self.pileup,
            sample_id=self.sample_id,
            chromosome=CHROM, sample_name=self.sample_name,
            sequencing_type=self.sequencing_type,
            output_base=self.output_directory + CHROM)
            for CHROM in CHROMs.iterkeys()]

if __name__ == "__main__":
    luigi.run()
