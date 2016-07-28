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
from parse_vcf_and_load import parse_vcf_and_load
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
    """Target describing the output of parsing a VCF/gVCF pair,
        and verification of the entries in the database
    """
    def __init__(self, vcf, chromosome, sample_id, output_base):
        self.vcf = vcf
        self.chromosome = chromosome
        self.sample_id = sample_id
        self.output_base = output_base
        self.novel_variants = output_base + ".novel_variants.txt"
        self.called_variants = output_base + ".calls.txt"
        self.variant_id_vcf = output_base + ".variant_id.vcf"
        self.matched_indels = output_base + ".matched_indels.txt"

    def exists(self):
        #print("checking files")
        for fn in (self.novel_variants, self.called_variants,
                   self.variant_id_vcf, self.matched_indels):
            if not os.path.isfile(fn):
                return False
        db = get_connection("dragen")
        try:
            vcf_lines = 0
            vcf_extra_calls = 0
            vcf_tabix = tabix.open(self.vcf)
            for vcf_fields in vcf_tabix.querys(self.chromosome):
                vcf_lines += 1
                vcf_extra_calls += vcf_fields[VCF_COLUMNS_DICT["ALT"]].count(",")
            vcf_variants_count = vcf_lines + vcf_extra_calls
            #print(vcf_variants_count)
            calls_line_count = get_num_lines_from_vcf(
                self.called_variants, header=False)
            #print(calls_line_count)
            if vcf_variants_count != calls_line_count:
                return False
            variant_id_vcf_line_count = get_num_lines_from_vcf(
                self.variant_id_vcf, header=False)
            #print(variant_id_vcf_line_count)
            if vcf_lines != variant_id_vcf_line_count:
                return False
            cur = db.cursor()
            cur.execute(GET_NUM_CALLS_FOR_SAMPLE.format(
                CHROM=self.chromosome, sample_id=self.sample_id))
            db_count = cur.fetchone()[0]
            #print(db_count)
            if db_count != vcf_variants_count:
                return False
            return True
        finally:
            if db.open:
                db.close()

    def open(self, mode="r"):
        return get_fh(self.called_variants, mode)

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

    def requires(self):
        return self.clone(CopyDataToScratch)

    def work(self):
        parse_vcf_and_load(
            vcf=self.copied_files_dict["vcf"],
            gvcf=self.copied_files_dict["gvcf"],
            pileup=self.copied_files_dict["pileup"],
            CHROM=self.chromosome, sample_id=self.sample_id,
            output_base=self.output_base, debug=self.debug)

    def output(self):
        return ParsedVCFTarget(
            vcf=self.copied_files_dict["vcf"], chromosome=self.chromosome,
            sample_id=self.sample_id, output_base=self.output_base)

class ImportSample(luigi.WrapperTask):
    """import a sample by requiring() a ParseVCF Task for each chromosome
    """
    vcf = luigi.InputFileParameter(description="the VCF to parse")
    gvcf = luigi.InputFileParameter(description="the gVCF to parse")
    pileup = luigi.InputFileParameter(description="the pileup to parse")
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
                self.sample_name, self.sequencing_type, self.capture_kit = row
            else:
                raise ValueError("sample_id {sample_id} does not exist".format(
                    sample_id=self.sample_id))
        finally:
            if db.open:
                db.close()
        self.output_directory = cfg.get("pipeline", "scratch_area").format(
            seqscratch=self.seqscratch, sample_name=self.sample_name,
            sequencing_type=self.sequencing_type.upper())

    def requires(self):
        return [self.clone(
            ParseVCF, chromosome=CHROM, sample_name=self.sample_name,
            sequencing_type=self.sequencing_type,
            output_base=self.output_directory + CHROM)
            for CHROM in CHROMs.iterkeys()]

if __name__ == "__main__":
    luigi.run()
