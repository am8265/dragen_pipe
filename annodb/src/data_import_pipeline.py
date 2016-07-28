#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Run pipeline to:
    1a. Annotate VCF
    1b. Split gVCF by chromosome
    2. Split VCF by chromosome
    3. Parse each chromosome's VCF & gVCF and import the data
    4. Archive appropriate data for sample
"""

import luigi
from luigi.contrib.sge import SGEJobTask
import os
import sys
import subprocess
import shlex
import operator
import re
import tabix
from shutil import copy
from collections import OrderedDict, Counter
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

class VCFTarget(luigi.Target):
    """Target that represents a VCF - checks for file existence and counts the
    number of lines
    """
    def __init__(self, vcf, vcf_lines):
        self.vcf = vcf
        self.vcf_lines = vcf_lines

    def exists(self):
        if not os.path.isfile(self.vcf):
            return False
        if self.vcf_lines != get_num_lines_from_vcf(self.vcf):
            return False
        return True

    def open(self, mode="r"):
        return get_fh(self.vcf, mode)

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

    def open(self):
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

class AnnotateVCF(SGEJobTask):
    """Annotate the specified VCF with SnpEff
    """
    vcf = luigi.InputFileParameter(description="the VCF to annotate")
    java = luigi.InputFileParameter(
        default=cfg.get("annotate", "java"),
        description="the path to the java executable for SnpEff")
    snpeff = luigi.InputFileParameter(
        default=cfg.get("annotate", "snpeff"),
        description="the path to the SnpEff jar file")
    genome_version = luigi.Parameter(
        default=cfg.get("annotate", "genome_version"),
        description="the Ensembl genome build to pass to SnpEff")
    snpeff_cfg = luigi.InputFileParameter(
        default=cfg.get("annotate", "snpeff_cfg"),
        description="the SnpEff config file to use")
    snpeff_options = luigi.Parameter(
        default=cfg.get("annotate", "snpeff_options"),
        description="additional parameters to pass to SnpEff")
    intervals = luigi.InputFileParameter(
        default=cfg.get("annotate", "intervals"),
        description="the intron/exon boundaries to annotate with")
    output_vcf = luigi.Parameter(
        default=None, description="the output VCF file path")
    poll_time = luigi.IntParameter(
        default=5, description="the time to wait to run qstat",
        significant=False)

    def __init__(self, *args, **kwargs):
        super(AnnotateVCF, self).__init__(*args, **kwargs)
        self.vcf_lines = get_num_lines_from_vcf(self.vcf)
        if self.output_vcf:
            if not os.path.isdir(os.path.dirname(self.output_vcf)):
                try:
                    os.makedirs(os.path.dirname(self.output_vcf))
                except:
                    raise OSError("couldn't create directory for output VCF")
        else:
            self.output_vcf = os.path.splitext(self.vcf)[0] + ".ann.vcf"

    def work(self):
        """perform the actual annotation of the VCF with SnpEff
        """
        cmd = ("{java} -Xmx5G -jar {snpeff} eff {genome_version} -c "
               "{snpeff_cfg} -interval {intervals} {snpeff_options} "
               "-o vcf {vcf}".format(**self.__dict__))
        with self.output().open("w") as vcf_out, \
                open(self.output_vcf + ".log", "w") as log_fh:
            p = subprocess.Popen(
                shlex.split(cmd), stdout=vcf_out, stderr=log_fh)
            p.wait()
        if p.returncode:
            raise subprocess.CalledProcessError(p.returncode, cmd)

    def output(self):
        return VCFTarget(vcf=self.output_vcf, vcf_lines=self.vcf_lines)

class SplitFileTarget(luigi.Target):
    """Target for testing existence and number of lines in each chromosome's
    file
    """
    def __init__(self, base_fn, split_files, header=True):
        self.base_fn = base_fn
        self.split_files = split_files
        self.header = header

    def exists(self):
        for fn in self.split_files.itervalues():
            if not os.path.isfile(fn):
                return False
        chromosome_lines = Counter()
        with get_fh(self.base_fn) as fh:
            if self.header:
                for line in fh:
                    if line.startswith("#CHROM"):
                        break
            for line in fh:
                CHROM = line.split("\t")[0]
                if CHROM in CHROMs:
                    chromosome_lines[CHROM] += 1
        print(chromosome_lines)
        for CHROM, line_count in chromosome_lines.iteritems():
            print(self.split_files[CHROM])
            if os.path.isfile(self.split_files[CHROM]):
                if (get_num_lines_from_vcf(
                    self.split_files[CHROM], header=False) != line_count):
                    return False
            else:
                return False
        return True

    def open(self, CHROM, mode="r"):
        return get_fh(self.split_files[CHROM], mode)

class SplitFileByChromosome(SGEJobTask):
    """Split a file by chromosome
    """
    fn = luigi.InputFileParameter(description="the file to parse")
    file_type = luigi.ChoiceParameter(
        choices=["vcf", "gvcf", "pileup"], description="the type of file")
    output_base = luigi.Parameter(
        default=None, description="specify a custom output file name")
    no_header = luigi.BoolParameter(
        default=False, description="specify if there is no header in the file")
    poll_time = luigi.IntParameter(
        default=5, description="the time to wait to run qstat",
        significant=False)

    def __init__(self, *args, **kwargs):
        super(SplitFileByChromosome, self).__init__(*args, **kwargs)
        if self.output_base:
            if re.search(r"\.g?vcf|\.pileup|\.bed", self.output_base):
                base_output = os.path.splitext(os.path.realpath(
                    self.output_base))[0]
            else:
                base_output = self.output_base
            if not os.path.isdir(os.path.dirname(base_output)):
                os.makedirs(os.path.dirname(base_output))
        else:
            base_output = os.path.splitext(self.fn)[0]

        suffix = "." + self.file_type
        self.split_files = OrderedDict([
            CHROM, base_output + "." + CHROM + suffix]
            for CHROM in CHROMs.iterkeys())

    def work(self):
        with get_fh(self.fn) as fh:
            last_CHROM = "-1"
            CHROM_idx = -1
            max_idx = len(CHROMs) - 1
            if not self.no_header:
                for line in fh:
                    if line.startswith("#CHROM"):
                        break
            for line in fh:
                CHROM = line.split("\t")[0]
                if CHROM != last_CHROM:
                    if last_CHROM != "-1":
                        out_fh.close()
                    if CHROM_idx == max_idx:
                        break
                    if CHROM in CHROMs:
                        new_CHROM_idx = CHROMs[CHROM]
                        if new_CHROM_idx < CHROM_idx:
                            raise ValueError(
                                "file is unsorted - chromosome {old} before "
                                "{new}".format(old=last_CHROM, new=CHROM))
                        else:
                            last_CHROM = CHROM
                            CHROM_idx = new_CHROM_idx
                            out_fh = self.output().open(CHROM, "w")
                if CHROM in CHROMs:
                    out_fh.write(line)
        if not out_fh.closed:
            out_fh.close()

    def output(self):
        return SplitFileTarget(
            base_fn=self.fn, split_files=self.split_files,
            header=not self.no_header)

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

    def work(self):
        parse_vcf_and_load(
            vcf=self.vcf, gvcf=self.gvcf, pileup=self.pileup,
            CHROM=self.chromosome, sample_id=self.sample_id,
            output_base=self.output_base)

    def output(self):
        return ParsedVCFTarget(
            vcf=self.vcf, chromosome=self.chromosome, sample_id=self.sample_id,
            output_base=self.output_base)

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
            ParseVCF, chromosome=CHROM, output_base=self.output_directory + CHROM)
            for CHROM in CHROMs.iterkeys()]

if __name__ == "__main__":
    luigi.run()
