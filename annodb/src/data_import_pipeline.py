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
from collections import OrderedDict, Counter
from parse_vcf_and_load import parse_vcf_and_load
from dragen_globals import *
from db_statements import *

cfg = get_cfg()
CHROMs = OrderedDict([CHROM, x] for x, CHROM in enumerate(
    cfg.get("ref", "CHROMs").split(",")))

def get_num_lines_from_vcf(vcf, header=True):
    """return the number of variant lines in the specified VCF, gzipped optional
    """
    count = 0
    with get_fh(vcf) as vcf_fh:
        if header:
            for line in vcf_fh:
                if line.startswith("#CHROM"):
                    break
        for _ in vcf_fh:
            count += 1
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

class SplitVCFTarget(luigi.Target):
    """Target for testing existence and number of lines in each chromosome's VCF
    """
    def __init__(self, base_vcf, split_vcfs):
        self.base_vcf = base_vcf
        self.split_vcfs = split_vcfs

    def exists(self):
        chromosome_lines = Counter()
        with get_fh(self.base_vcf) as vcf_fh:
            for line in vcf_fh:
                if line.startswith("#CHROM"):
                    break
            for line in vcf_fh:
                CHROM = line.split("\t")[0]
                if CHROM in CHROMs:
                    chromosome_lines[CHROM] += 1
        for CHROM, line_count in chromosome_lines.iteritems():
            if os.path.isfile(self.split_vcfs[CHROM]):
                if (get_num_lines_from_vcf(
                    self.split_vcfs[CHROM], header=False) != line_count):
                    return False
            else:
                return False
        return True

    def open(self, CHROM, mode="r"):
        return get_fh(self.split_vcfs[CHROM], mode)

class SplitVCFByChromosome(SGEJobTask):
    """Split a VCF by chromosome
    """
    vcf = luigi.InputFileParameter(description="the VCF to parse")
    gvcf = luigi.BoolParameter(
        default=False, description="whether the VCF is a gVCF or not")
    output_base = luigi.Parameter(
        default=None, description="specify a custom output file name")
    poll_time = luigi.IntParameter(
        default=5, description="the time to wait to run qstat",
        significant=False)

    def __init__(self, *args, **kwargs):
        super(SplitVCFByChromosome, self).__init__(*args, **kwargs)
        if self.output_base:
            if re.search(r"\.g?vcf", self.output_base):
                base_output = os.path.splitext(os.path.relapath(
                    self.output_base))[0]
            else:
                base_output = self.output_base
            if not os.path.isdir(os.path.dirname(base_output)):
                os.makedirs(os.path.dirname(base_output))
        else:
            base_output = os.path.splitext(self.vcf)[0]

        suffix = ".gvcf" if self.gvcf else ".vcf"
        self.split_vcfs = OrderedDict([
            CHROM, base_output + "." + CHROM + suffix]
            for CHROM in CHROMs.iterkeys())

    def work(self):
        with open(self.vcf) as vcf_fh:
            last_CHROM = "-1"
            CHROM_idx = -1
            max_idx = len(CHROMs) - 1
            for line in vcf_fh:
                if line.startswith("#CHROM"):
                    break
            for line in vcf_fh:
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
                                "VCF is unsorted - chromosome {old} before {new}".
                                format(old=last_CHROM, new=CHROM))
                        else:
                            last_CHROM = CHROM
                            CHROM_idx = new_CHROM_idx
                            out_fh = self.output().open(CHROM, "w")
                if CHROM in CHROMs:
                    out_fh.write(line)
        if not out_fh.closed:
            out_fh.close()

    def output(self):
        return SplitVCFTarget(base_vcf=self.vcf, split_vcfs=self.split_vcfs)

class ParsedVCFTarget(luigi.Target):
    """Target describing the output of parsing a VCF/gVCF pair,
        and verification of the entries in the database
    """
    def __init__(self, vcf, chromosome, sample_id):
        self.vcf = vcf
        self.chromosome = chromosome
        self.sample_id = sample_id
        self.called_variants = self.vcf + ".calls.txt"

    def exists(self):
        if not os.path.isfile(self.called_variants):
            return False
        db = get_connection("dragen")
        try:
            vcf_line_count = get_num_lines_from_vcf(self.vcf, header=False)
            calls_line_count = get_num_lines_from_vcf(
                self.called_variants, header=False)
            if vcf_line_count != calls_line_count:
                return False
            cur = db.cursor()
            cur.execute("SELECT COUNT(variant_id) FROM called_variant_chr{CHROM} "
                        "WHERE sample_id = {sample_id}".format(
                            CHROM=self.chromosome, sample_id=self.sample_id))
            db_count = cur.fetchone()[0]
            if db_count != vcf_line_count:
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
    chromosome = luigi.Parameter(
        description="the chromosome that's being processed")
    sample_id = luigi.NumericalParameter(
        left_op=operator.le, description="the sample_id for this sample")

    def __init__(self, *args, **kwargs):
        super(SGEJobTask, self).__init__(*args, **kwargs)
        if self.chromosome not in CHROMs:
            raise ValueError("invalid chromosome specified: {chromosome}".format(
                chromosome=self.chromosome))

    def requires(self):
        return self.clone(SplitVCFByChromosome)

    def work(self):
        parse_vcf_and_load(
            vcf=self.vcf, gvcf=self.gvcf, chromosome=self.chromosome,
            sample_id=self.sample_id)

    def output(self):
        return ParsedVCFTarget(
            vcf=self.vcf, chromosome=self.chromsome, sample_id=self.sample_id)

if __name__ == "__main__":
    luigi.run()
