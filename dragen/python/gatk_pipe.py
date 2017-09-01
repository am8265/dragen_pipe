#!/nfs/goldstein/software/python2.7.7/bin/python2.7

import os
import sys
import subprocess
import luigi
import MySQLdb
import warnings
import tabix
import grp
import tarfile
from random import random
from time import sleep
from glob import glob
from numpy import isclose
from shutil import rmtree
from collections import Counter
from getpass import getuser
from pwd import getpwuid
from gzip import open as gopen
from dragen_db_statements import *
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.realpath(__file__)))), "import", "src"))
from waldb_globals import *
from check_avere import check_avere

def owner(fn):
    """Return the user name of the owner of the specified file
    """
    if os.path.exists(fn):
        return getpwuid(os.stat(fn).st_uid).pw_name
    else:
        raise OSError("{fn} does not exist!".format(fn=fn))

def getDBIDMaxPrepID(pseudo_prepid):
    """ Function to get max prepid for a given pseudo_prepid
    from the database table pseudo_prepid

    pseudo_prepid : str : the pseudo_prepid value """

    db = get_connection("seqdb")
    try:
        cur = db.cursor()
        cur.execute(GET_DBID_MAX_PREPID.format(
            pseudo_prepid=pseudo_prepid))
        results = cur.fetchone()
        dbid,prepid = results
    finally:
        if db.open:
            db.close()

    return dbid,prepid

def getCaptureKit(pseudo_prepid):
    """ Get the capturekit for a given pseudo_prepid

    pseudo_prepid : str : the pseudo_prepid value """

    db = get_connection("seqdb")
    try:
        cur = db.cursor()
        cur.execute(GET_CAPTURE_KIT.format(
            pseudo_prepid=pseudo_prepid))
        capture_kit = cur.fetchone()
    finally:
        if db.open:
            db.close()
    return capture_kit

def getPanelOfNormals(captureKit,sample_type):
    """ Get the panel of normals for CNV normalization

    captureKit : str ; the capturekit , e.g. Roche
    sample_type : str ; the sequencing type, e.g exome """

    db = get_connection("seqdb")
    try:
        cur = db.cursor()
        cur.execute(GET_PANEL_OF_NORMALS.format(
            captureKit=captureKit,
            sample_type=sample_type))
        panelOfNormals = cur.fetchone()
    finally:
        if db.open:
            db.close()
    return panelOfNormals

class programs(luigi.Config):
    """ Configuration class for instantiating parameters for this pipeline
    the values are read from luigi.cfg in the current folder """
    bedtools = luigi.InputFileParameter()
    bgzip = luigi.InputFileParameter()
    gatk = luigi.InputFileParameter()
    gatk4 = luigi.InputFileParameter()
    java = luigi.InputFileParameter()
    picard = luigi.InputFileParameter()
    pypy = luigi.InputFileParameter()
    snpEff = luigi.InputFileParameter()
    snpSift = luigi.InputFileParameter()
    tabix = luigi.InputFileParameter()
    verifybamid = luigi.InputFileParameter()
    vcffilter = luigi.InputFileParameter()
    clineff = luigi.InputFileParameter()
    python_path = luigi.InputFileParameter()

class locations(luigi.Config):
    base = luigi.InputDirectoryParameter()

class pipeline_files(luigi.Config):
    ccds_bed_file = luigi.InputFileParameter()
    coverage_binner = luigi.InputFileParameter()
    gq_binner = luigi.InputFileParameter()
    snpEff_cfg = luigi.InputFileParameter()
    clineff_cfg = luigi.InputFileParameter()
    target_file = luigi.InputFileParameter()
    target_file_X = luigi.InputFileParameter()
    target_file_Y = luigi.InputFileParameter()
    transpose_awk = luigi.InputFileParameter()
    subset_vcf_awk = luigi.InputFileParameter()

class gatk_resources(luigi.Config):
    contam1000g_vcf = luigi.InputFileParameter()
    ref_genome = luigi.InputFileParameter()
    seqdict_file = luigi.InputFileParameter()
    interval = luigi.InputFileParameter()
    g1000 = luigi.InputFileParameter()
    hapmap = luigi.InputFileParameter()
    Mills1000g = luigi.InputFileParameter()
    omni = luigi.InputFileParameter()
    dbSNP = luigi.InputFileParameter()
    annotatedbSNP = luigi.InputFileParameter()
    exac = luigi.InputFileParameter()

class variables(luigi.Config):
    create_targetfile = luigi.BooleanParameter()
    max_mem = luigi.IntParameter()
    cnv_target_padding = luigi.IntParameter(
        description="Num. of bp to pad each target during "
        "GATK's CalculateTargetCoverage CNV task")
    block_size = luigi.IntParameter(description='The block size over which to do the binning')

class qc_metrics(luigi.Config):
    """ Parse in database field names for qc metrics
    from the qc_metrics section of the config file """

    total_reads = luigi.Parameter()
    pct_reads_aligned = luigi.Parameter()
    pct_mismatch_rate = luigi.Parameter()
    capture_specificity = luigi.Parameter()
    mean_cvg = luigi.Parameter()
    median_cvg = luigi.Parameter()
    mean_median_cvg = luigi.Parameter()
    pct_bases5X = luigi.Parameter()
    pct_bases10X = luigi.Parameter()
    pct_bases15X = luigi.Parameter()
    pct_bases20X = luigi.Parameter()
    mean_ccds_cvg = luigi.Parameter()
    pct_ccds_bases5X = luigi.Parameter()
    pct_ccds_bases10X = luigi.Parameter()
    pct_ccds_bases15X = luigi.Parameter()
    pct_ccds_bases20X = luigi.Parameter()
    mean_X_cvg = luigi.Parameter()
    mean_Y_cvg = luigi.Parameter()
    xy_avg_ratio = luigi.Parameter()
    pct_duplicate_reads = luigi.Parameter()
    total_snps = luigi.Parameter()
    pct_dbsnp_snps = luigi.Parameter()
    total_indels = luigi.Parameter()
    pct_dbsnp_indels = luigi.Parameter()
    titv = luigi.Parameter()
    homhet_ratio = luigi.Parameter()
    snv_homhet_ratio = luigi.Parameter()
    indel_homhet_ratio = luigi.Parameter()
    x_homhet_ratio = luigi.Parameter()
    contamination_value = luigi.Parameter()
    concordance = luigi.Parameter()

# create one dict that contains all these parameters for ease of use

class GATKFPipelineTask(GATKPipelineTask):
    task_name_format = "{task_family}.{sample_name}.{pseudo_prepid}"
    poll_time = 120
    def __init__(self, *args, **kwargs):
        self.config_parameters = {}
        for cls in (programs, locations, pipeline_files, gatk_resources, variables,
                    qc_metrics):
            self.config_parameters.update(cls().__dict__)
        super(GATKFPipelineTask, self).__init__(*args, **kwargs)

    def format_string(self, s, **kwargs):
        """Format the string s with any parameters in config_parameters, self, and kwargs
        (priority if overlap is given to kwargs then self)
        """
        self.config_parameters.update(self.__dict__)
        self.config_parameters.update(kwargs)
        return s.format(**self.config_parameters)

class RealignerTargetCreator(GATKFPipelineTask):
    priority = 1 # run before other competing steps with no requirement
    n_cpu = 4
    def pre_shell_commands(self):
        self.commands = [self.format_string(
            "{java} -Xmx{max_mem}g -jar {gatk} -R {ref_genome} -T RealignerTargetCreator "
            "-I {scratch_bam} -o {interval_list} -known {Mills1000g} "
            "-known {dbSNP} --log_to_file {log_file} -nt {n_cpu}")]

class IndelRealigner(GATKFPipelineTask):
    def pre_shell_commands(self):
        self.commands = [self.format_string(
            "{java} -Xmx{max_mem}g -jar {gatk} -R {ref_genome} -T IndelRealigner "
            "-I {scratch_bam} -o {realn_bam} "
            "--log_to_file {log_file} -targetIntervals {interval_list} "
            "-maxReads 10000000 -maxInMemory 450000 -known {Mills1000g} "
            "-known {dbSNP}")]

    def _run_post_success(self):
        os.remove(self.scratch_bam)

    def requires(self):
        return self.clone(RealignerTargetCreator)

class BaseRecalibrator(GATKFPipelineTask):
    n_cpu = 4
    def pre_shell_commands(self):
        self.commands = [self.format_string(
            "{java} -Xmx{max_mem}g -jar {gatk} -R {ref_genome} "
            "-T BaseRecalibrator -I {realn_bam} "
            "-o {recal_table} --log_to_file {log_file} "
            "-knownSites {Mills1000g} -knownSites {dbSNP} -nct {n_cpu}")]

    def requires(self):
        return self.clone(IndelRealigner)

class PrintReads(GATKFPipelineTask):
    n_cpu = 4
    def pre_shell_commands(self):
        # --disable_indel_quals are necessary to remove BI and BD tags in the bam file
        self.commands = [self.format_string(
            "{java} -Xmx{max_mem}g -jar {gatk} -R {ref_genome} -T PrintReads "
            "-I {realn_bam} --log_to_file {log_file} "
            "--disable_indel_quals -BQSR {recal_table} -o {recal_bam} "
            "-nct {n_cpu}")]

    def _run_pos_success(self):
        os.remove(self.realn_bam)

    def requires(self):
        return self.clone(BaseRecalibrator)

class HaplotypeCaller(GATKFPipelineTask):
    priority = 1 # run before other steps needing the recalibrated BAM
    n_cpu = 6
    def pre_shell_commands(self):
        self.commands = [self.format_string(
            "{java} -Xmx{max_mem}g -jar {gatk} -R {ref_genome} -T HaplotypeCaller "
            "-L {interval} -I {recal_bam} -o {gvcf} --log_to_file {log_file} "
            "-stand_call_conf 20 -stand_emit_conf 20 --emitRefConfidence GVCF "
            "-GQB 5 -GQB 15 -GQB 20 -GQB 60 --variant_index_type LINEAR "
            "--variant_index_parameter 128000 --dbsnp {dbSNP} -nct {n_cpu}")]

    def requires(self):
        return self.clone(PrintReads)

class GenotypeGVCFs(GATKFPipelineTask):
    priority = 1
    def pre_shell_commands(self):
        self.commands = [self.format_string(
            "{java} -Xmx{max_mem}g -jar {gatk} -R {ref_genome} -T GenotypeGVCFs "
            "-L {interval} -o {vcf} --log_to_file {log_file} "
            "-stand_call_conf 20 -stand_emit_conf 20 -V {gvcf}")]

    def requires(self):
        return self.clone(HaplotypeCaller)

class SelectVariantsSNP(GATKFPipelineTask):
    priority = 1
    def pre_shell_commands(self):
        self.commands = [self.format_string(
            "{java} -Xmx{max_mem}g -jar {gatk} -R {ref_genome} -T SelectVariants "
            "-L {interval} -V {vcf} --log_to_file {log_file} -selectType SNP "
            "-o {snp_vcf}")]

    def requires(self):
        return self.clone(GenotypeGVCFs)

class SelectVariantsINDEL(GATKFPipelineTask):
    def pre_shell_commands(self):
        self.commands = [self.format_string(
            "{java} -Xmx{max_mem}g -jar {gatk} -R {ref_genome} -T SelectVariants "
            "-L {interval} -V {vcf} --log_to_file {log_file} -selectType INDEL "
            "-o {indel_vcf}")]

    def requires(self):
        return self.clone(GenotypeGVCFs)

class VariantRecalibratorSNP(GATKFPipelineTask):
    def pre_shell_commands(self):
        ## Running both Exomes and Genomes the same way, i.e. we will exlcude DP, omni is set to be a truth set, also added in ExAC SNPs as a training
        ## set which can contain both TP and FP with the same prior likelihood as the 1000G training set. 
        ## See these links : 
        ## For exomes : 
        ## https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php
        ## For genomes : 
        ## https://software.broadinstitute.org/gatk/guide/article?id=1259
        self.commands = [self.format_string(
            "{java} -Xmx{max_mem}g -jar {gatk} -R {ref_genome} -T VariantRecalibrator "
            "-L {interval} --log_to_file {log_file} --input {snp_vcf} "
            "-an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum "
            "-mode SNP --maxGaussians 4 -tranche 100.0 -tranche 99.9 "
            "-tranche 99.0 -tranche 90.0 -recalFile {snp_recal} "
            "-tranchesFile {snp_tranches} -rscriptFile {snp_rscript} "
            "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} "
            "-resource:omni,known=false,training=true,truth=true,prior=12.0 {omni} "
            "-resource:1000G,known=false,training=true,truth=false,prior=10.0 {g1000} "
            "-resource:ExAc,known=false,training=true,truth=false,prior=10.0 {exac} "
            "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbSNP}")]

    def requires(self):
        return self.clone(SelectVariantsSNP)

class VariantFiltrationSNP(GATKFPipelineTask):
    def pre_shell_commands(self):
        self.commands = [self.format_string(
            '{java} -Xmx{max_mem}g -jar {gatk} -R {ref_genome} -T VariantFiltration '
            '-L {interval} -V {snp_vcf} --filterExpression '
            '"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || '
            'ReadPosRankSum < -8.0" --filterName "SNP_filter" '
            '--log_to_file {log_file} -o {snp_filtered}')]

    def requires(self):
        return self.clone(SelectVariantsSNP)

class VariantRecalibratorINDEL(GATKFPipelineTask):
    def pre_shell_commands(self):
        self.commands = [self.format_string(
            "{java} -Xmx{max_mem}g -jar {gatk} -R {ref_genome} "
            "-T VariantRecalibrator -L {interval} --log_to_file {log_file} "
            "--input {indel_vcf} -an QD -an FS -an SOR -an MQRankSum "
            "-an ReadPosRankSum -mode INDEL --maxGaussians 4 "
            "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 "
            "-recalFile {indel_recal} -tranchesFile {indel_tranches} "
            "-resource:mills,known=true,training=true,truth=true,prior=12.0 {Mills1000g}")]

    def requires(self):
        return self.clone(SelectVariantsINDEL)

class ApplyRecalibrationSNP(GATKFPipelineTask):
    def pre_shell_commands(self):
        commands = [self.format_string(
            "{java} -Xmx{max_mem}g -jar {gatk} -R {ref_genome} -T ApplyRecalibration "
            "-L {interval} --log_to_file {log_file} -input {snp_vcf} "
            "-tranchesFile {snp_tranches} -recalFile {snp_recal} "
            "-o {snp_filtered} --ts_filter_level 90.0 -mode SNP")]

    def requires(self):
        return self.clone(VariantRecalibratorSNP)

class ApplyRecalibrationINDEL(GATKFPipelineTask):
    def pre_shell_commands(self):
        self.commands = [self.format_string(
            "{java} -Xmx{max_mem}g -jar {gatk} -R {ref_genome} -T ApplyRecalibration "
            "-L {interval} --log_to_file {log_file} -input {indel_vcf} "
            "-tranchesFile {indel_tranches} -recalFile {indel_recal} "
            "-o {indel_filtered} --ts_filter_level 90.0 -mode INDEL")]

    def requires(self):
      return self.clone(VariantRecalibratorINDEL)

class VariantFiltrationINDEL(GATKFPipelineTask):
    def pre_shell_commands(self):
        self.commands = [self.format_string(
            '{java} -Xmx{max_mem}g -jar {gatk} -R {ref_genome} -T VariantFiltration '
            '-L {interval} -V {indel_vcf} --filterExpression '
            '"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName '
            '"INDEL_filter" --log_to_file {log_file} -o {indel_filtered}')]

    def requires(self):
      return self.clone(SelectVariantsINDEL)

class CombineVariants(GATKFPipelineTask):
    def pre_shell_commands(self):
        """Merges SNP and INDEL vcfs.  Using the filtered SNP vcf header as a
        base, variant type/sample type specific ##FILTERs are added to the header.
        After finishing reading through the header the SNP vcf is read again with 
        only variants being outputted.  The same happens with the INDEL vcf.
        The resulting vcf is then sorted producing the analysisReady.vcf.
        """
        filter_flag = 0
        info_flag = 0
        
        with open(self.tmp_vcf, "w") as vcf_out, open(self.snp_filtered) as header:
            for line in header.readlines():
                if line.startswith("#"):
                    if line.startswith("##FILTER") and filter_flag == 0 :
                        filter_flag = 1
                        #SNP specific FILTER whether using VQSR or snp filtering
                        if self.sample_type in ("EXOME", "GENOME"):            ## Will a 90.0to99.00 tranche be needed in the header ?
                            vcf_out.write('##FILTER=<ID=VQSRTrancheSNP90.00to99.00,Description="Truth sensitivity tranche level for SNP model"\n')
                            vcf_out.write('##FILTER=<ID=VQSRTrancheSNP99.00to99.90,Description="Truth sensitivity tranche level for SNP model"\n')
                            vcf_out.write('##FILTER=<ID=VQSRTrancheSNP99.90to100.00+,Description="Truth sensitivity tranche level for SNP model"\n')
                            vcf_out.write('##FILTER=<ID=VQSRTrancheSNP99.90to100.00,Description="Truth sensitivity tranche level for SNP model"\n')
                        else:
                            vcf_out.write('##FILTER=<ID=SNP_filter,Description="QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0">\n')

                        #Indel specific filters whether using VQSR or indel filtering
                        if self.sample_type == "GENOME":
                            vcf_out.write('##FILTER=<ID=VQSRTrancheINDEL90.00to99.00,Description="Truth sensitivity tranche level for INDEL model"\n')
                            vcf_out.write('##FILTER=<ID=VQSRTrancheINDEL99.00to99.90,Description="Truth sensitivity tranche level for INDEL model"\n')
                            vcf_out.write('##FILTER=<ID=VQSRTrancheINDEL99.90to100.00+,Description="Truth sensitivity tranche level for INDEL model"\n')
                            vcf_out.write('##FILTER=<ID=VQSRTrancheINDEL99.90to100.00,Description="Truth sensitivity tranche level for INDEL model"\n')
                        else:
                            vcf_out.write('##FILTER=<ID=INDEL_filter,Description="QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0">\n')

                    if line.startswith("##INFO") and info_flag == 0:
                        info_flag = 1
                        ## Adding this wierd SAO info tag since there are rare cases where a variant has this information, but it is absent in the header
                        ## causing read backed phasing to crash.
                        vcf_out.write("""##INFO=<ID=SAO,Number=1,Type=Integer,Description="Variant Allele Origin: 0. unspecified, 1. Germline, 2. Somatic, 3. Both">\n""")
                    vcf_out.write(line)
                else:
                    break

            with open(self.snp_filtered) as snps:
                for snp in snps.readlines():
                    if not snp.startswith("#"):
                        vcf_out.write(snp)
            with open(self.indel_filtered) as indels:
                for indel in indels.readlines():
                    if not indel.startswith("#"):
                        vcf_out.write(indel)

        self.commands.append(self.format_string(
            "{java} -jar {picard} SortVcf I={tmp_vcf} O={final_vcf}"))
        self.commands.append(self.format_string("{bgzip} -f {final_vcf}"))
        self.commands.append(self.format_string(
            "{tabix} -f {final_vcf_gz}"))

    def requires(self):
        if self.sample_type in ("EXOME", "GENOME"):
            yield self.clone(ApplyRecalibrationSNP)
            if self.sample_type == "EXOME":
                yield self.clone(VariantFiltrationINDEL)
            else:
                yield self.clone(ApplyRecalibrationINDEL)
        elif self.sample_type == "CUSTOM_CAPTURE":
            yield self.clone(VariantFiltrationSNP)
            yield self.clone(VariantFiltrationINDEL)
        else:
            raise ValueError(
                "Sample type: {} not supported".format(self.sample_type))

class RBP(GATKFPipelineTask):
    priority = 1
    def pre_shell_commands(self):
        self.commands = [self.format_string(
            "{java} -jar {gatk} -T ReadBackedPhasing -R {ref_genome} "
            "-I {recal_bam} --variant {final_vcf_gz} -o {phased_vcf} "
            "--phaseQualityThresh 20.0 --maxGenomicDistanceForMNP 2 "
            "--enableMergePhasedSegregatingPolymorphismsToMNP "
            "-U ALLOW_SEQ_DICT_INCOMPATIBILITY")]

    def requires(self):
        return self.clone(CombineVariants), self.clone(PrintReads)

class FixMergedMNPInfo(GATKFPipelineTask):
    """ The MNPs phased by RBP are missing the DP,AD and other annotation info
    this task which check variants which are missing these information, fetch
    the value for corresponding to that variant site from the vcf prior to RBP
    and add these information to the fixed vcf """
    def post_shell_commands(self):
        self.fix_phased_vcf()
        os.remove(self.phased_vcf)

    def fix_phased_vcf(self):
        """ Fix the missing DP and AD fields for the phased variants
        be gzipped and tabix indexed
        """
        self.vcf_tabix = tabix.open(self.final_vcf_gz)
        with open(self.phased_vcf, "r") as PHASED_VCF, open(self.fixed_vcf, "w") as OUT:
            for line in PHASED_VCF:
                line = line.strip("\n")
                fields = line.split("\t")
                if line.startswith("#"):
                    print >> OUT, line
                elif "DP" in fields[VCF_COLUMNS_DICT["FORMAT"]].split(":"):
                    print >> OUT, line
                else:
                    lREF = len(fields[VCF_COLUMNS_DICT["REF"]])
                    lALT = len(fields[VCF_COLUMNS_DICT["ALT"]])
                    if (lREF > 1 and lREF == lALT and "," not in 
                        fields[VCF_COLUMNS_DICT["ALT"]]):
                        # valid MNP
                        # get DP,AD_ref,AD_alt
                        chrom = fields[VCF_COLUMNS_DICT["CHROM"]]
                        pos = int(fields[VCF_COLUMNS_DICT["POS"]])
                        ad, dp = self.get_ad_dp_from_vcf(chrom, pos)
                        annotations = self.get_gatk_annotations_from_vcf(chrom, pos)
                        # Construct a new vcf line
                        new_line = self.construct_fixed_vcf_line(line,[ad,dp,annotations])
                        print >> OUT,new_line
                    else: # A non MNP with no DP, raise an exception
                        raise ValueError("A non MNP encountered with no DP in info field")

    def construct_fixed_vcf_line(self, line, missing_vals):
        """ Construct the fixed line 

        line : str ; the vcf line
        missing_vals : [str,str,str] ; the missing dp,ad and
        gatk annotations like ExcessHet,FS,VQSLOD,etc.

        returns : the line with the missing values fixed
        """

        info_fields = line.split("\t")[-2].split(":")
        vals = line.split("\t")[-1].split(":")

        ## Insert the new tags at the correct position
        info_fields.insert(1, "AD")
        info_fields.insert(2, "DP")

        vals.insert(1,missing_vals[0])
        vals.insert(2,missing_vals[1])

        new_info = ":".join(info_fields)
        new_vals = ":".join(vals)

        ## Fix the missing GATK annotations
        annotations = line.split("\t")[7]
        new_annotations = annotations + ";" + missing_vals[2]

        contents = line.split("\t")
        contents[-1] = new_vals
        contents[-2] = new_info
        contents[7] = new_annotations

        ## Return the fixed line
        return "\t".join(contents)

    def get_gatk_annotations_from_vcf(self, chrom, pos):
        """ Get GATK specific annotations like FS,VQSLOD,ExcessHet
        which are absent in MNPs after RBP has been run

        chrom : str ; the chromosome number
        pos : int ; the position for the variant

        returns : The annotations with the AC,AF,AN removed
        """
        return (";".join(list(self.vcf_tabix.query(chrom,pos,pos))[0][7].split(";")[3:]))

    def get_correct_record(self,records,pos):
        """ Return the record having its position equal to pos

        record : a list of lists ; obtained using a tabix tabix query
        pos : int ; the exact variant position we want

        returns : list ; the record corresponding to the exact pos """

        for record in records:
            if int(record[1]) == pos:
                return record

        ## The position was not found
        message="The variant site {0} could not be found in the vcf prior to RBP".format(pos)
        raise Exception(message)

    def get_ad_dp_from_vcf(self, chrom, pos):
        """ Get DP, AD_REF, AD_ALT
        from a VCF file

        chrom : str ; the chromosome number
        pos : int ; the position for the variant

        returns : [dp,ad_ref,ad_alt]
        """
        records = list(self.vcf_tabix.query(chrom, pos, pos))
        # The above query returns all variant sites
        # spanning the pos which we have queried for
        # We only want the exact pos site, so check
        # for that condition
        record = self.get_correct_record(records,pos)
        info_fields = record[-2].split(":")
        values = record[-1].split(":")
        found = 0
        for i in range(0,len(info_fields)):
            if info_fields[i] == "DP":
                dp = values[i]
                found+=1
            if info_fields[i] == "AD":
                ad = values[i]  # How to deal with multi allelic sites ? 
                found+=1

        if found != 2:
            raise Exception("Could not parse the vcf to get DP/AD info for "
                            "the phased variants")

        return [ad,dp]

    def requires(self):
        return self.clone(RBP)

class AnnotateVCF(GATKFPipelineTask):
    def pre_shell_commands(self):
        self.shell_options.update(
            {"stdout":None, "stderr":self.log_file, "shell":True})
        self.commands = [self.format_string(
            "/nfs/goldstein/software/jdk1.8.0_05/bin/java -Xmx6g -jar {clineff} "
            "-c {clineff_cfg} -v -db {annotatedbSNP} GRCh37.87 {fixed_vcf} "
            "> {annotated_vcf}"),
            self.format_string("{bgzip} -f {annotated_vcf}"),
            self.format_string("{tabix} -f {annotated_vcf_gz}")]

    def requires(self):
        return self.clone(FixMergedMNPInfo)

class SubsetVCF(GATKFPipelineTask):
    """ For SRR/Any future deemed to be low quality sequenced samples, subset the vcf to only the capturekit regions """
    def pre_shell_commands(self):
        self.shell_options.update({"stdout":None, "stderr":self.log_file})
        self.original_vcf_gz = "{0}/{1}.{2}.analysisReady.annotated.original.vcf.gz".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.restricted_vcf = "{0}/{1}.{2}.analysisReady.annotated.restricted.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.restricted_vcf_gz = "{0}/{1}.{2}.analysisReady.annotated.restricted.vcf.gz".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.subset_cmd = "{0} {1} -B {2} -h > {3}".format(
            self.config_parameters["tabix"],self.annotated_vcf_gz,self.capture_kit_bed,
            self.restricted_vcf)
        self.bgzip_cmd = "{0} -f {1}".format(self.config_parameters["bgzip"],self.restricted_vcf)
        self.rename1 = "mv {0} {1}".format(self.annotated_vcf_gz,self.original_vcf_gz)
        self.rename2 = "mv {0} {1}".format(self.restricted_vcf_gz,self.annotated_vcf_gz)
        self.tabix_cmd1 = "{0} -f {1}".format(self.config_parameters["tabix"],self.original_vcf_gz)
        self.tabix_cmd2 = "{0} -f {1}".format(self.config_parameters["tabix"],self.annotated_vcf_gz)
        self.shell_options.update({"stdout":None, "shell":True})
        self.commands = [
            self.subset_cmd, self.bgzip_cmd, self.rename1, self.rename2,
            self.tabix_cmd1, self.tabix_cmd2]

    def requires(self):
        return self.clone(AnnotateVCF)

class ArchiveSample(GATKFPipelineTask):
    """ Archive samples on Amplidata """
    n_cpu = 7
    def pre_shell_commands(self):
        # wait for Avere to cooperate before loading
        while True:
            if check_avere():
                break
            else:
                sleep(25 + 10 * random())
        if (self.sample_name.upper().startswith('PGMCLIN') or
            self.sample_name.upper().startswith('PGMVIP')):
            self.base_dir = os.path.join(
                "/nfs/pgmclin/ALIGNMENT/BUILD37/DRAGEN/",
                self.sample_type, self.name_prep)
        else:
            self.base_dir = os.path.join(
                self.config_parameters["base"], self.sample_type, self.name_prep)
        self.script_dir = os.path.join(self.scratch_dir, "scripts")
        self.pipeline_tarball = (
            "{scratch_dir}/{name_prep}.pipeline_data.tar.gz".format(
                scratch_dir=self.scratch_dir, name_prep=self.name_prep))
        self.cvg_tarball = (
            "{cov_dir}/coverage.tar.gz".format(cov_dir=self.cov_dir))
        self.gq_tarball =  (
            "{gq_dir}/gq.tar.gz".format(gq_dir=self.gq_dir))
        if not os.path.isdir(self.base_dir):
            os.makedirs(self.base_dir)
        with tarfile.open(self.pipeline_tarball, "w:gz") as tar:
            for d in (self.script_dir, self.log_dir):
                tar.add(d, arcname=os.path.basename(d))
            for txt in glob("{scratch_dir}/*.txt".format(scratch_dir=self.scratch_dir)):
                tar.add(txt, arcname=os.path.basename(txt))
        with tarfile.open(self.cvg_tarball, "w:gz") as tar:
            for cvg_file in glob(
                "{cov_dir}/*.txt".format(cov_dir=self.cov_dir)):
                tar.add(cvg_file, arcname=os.path.basename(cvg_file))
        with tarfile.open(self.gq_tarball, "w:gz") as tar:
            for gq_file in glob(
                "{gq_dir}/*.txt".format(gq_dir=self.gq_dir)):
                tar.add(gq_file, arcname=os.path.basename(gq_file))

        cmd = ["rsync", "-grlt", "--inplace", "--partial", "{pipeline_tarball}",
               "{recal_bam}", "{recal_bam_index}", "{annotated_vcf_gz}",
               "{annotated_vcf_gz_index}", "{gvcf}", "{gvcf_index}",
               "{cvg_tarball}", "{gq_tarball}"]
        if self.sample_type == "GENOME":
            if os.path.isfile(self.original_vcf_gz):
                cmd.append("{original_vcf_gz}")
        cmd.append("{base_dir}")
        self.commands = [self.format_string(" ".join(cmd))]

    def post_shell_commands(self):
        ## Change folder permissions of base directory
        uid = os.getuid()
        gid = grp.getgrnam("bioinfo").gr_gid
        os.chown(self.base_dir, uid, gid)
        user = getuser()
        for root, dirs, files in os.walk(self.base_dir):
            for d in dirs:
                d = os.path.join(root, d)
                if not os.path.islink(d) and owner(d) == user:
                    os.chmod(d, 0775)
                    os.chown(d, uid, gid)
            for f in files:
                f = os.path.join(root, f)
                if not os.path.islink(f) and owner(f) == user:
                    os.chmod(f, 0664)
                    os.chown(f, uid, gid)
        if os.path.isdir(self.scratch_dir):
            rmtree(self.scratch_dir)

        ## Update the AlignSeqFile loc to the final archive location
        location = "{0}/{1}".format(
            self.config_parameters["base"], self.sample_type)
        db = get_connection("seqdb")
        try:
            cur = db.cursor()
            update_alignseqfile(cur, location, self.pseudo_prepid)
            db.commit()
        finally:
            if db.open:
                db.close()

    def requires(self):
        if self.sample_name.upper().startswith('SRR'):
            yield self.clone(SubsetVCF)
        yield self.clone(UpdateSeqdbMetrics)

def get_productionvcf(pseudo_prepid, sample_name, sample_type):
    """
    Query seqdbClone to get AlignSeqFileLoc
    Extend that to search for the production vcf
    file, tabix index it, if not already.

    Args :    pseudo_prepid; Int; An identifier for the sample in the database
              sample_name : String; CHGVID, a string identifier for the 
                            sequencing sample
              sample_type : String; The sequencing type i.e. Exome,Genomee,etc.

    Returns : String; Path to the production vcf file
              None if vcf was not found or if the sample
              was not run through the production pipeline
              before
    """

    DBID,prepID = getDBIDMaxPrepID(pseudo_prepid)
    query_statement = (
        """ SELECT AlignSeqFileLoc FROM seqdbClone WHERE"""
        """ CHGVID = '{0}' AND seqType = '{1}' AND """
        """ pseudo_prepid = '{2}'""".format(
            sample_name, sample_type, pseudo_prepid))
    db = get_connection("seqdb")
    try:
        cur = db.cursor()
        cur.execute(query_statement)
        db_val = cur.fetchall()
        if len(db_val) == 0: # If the query returned no results
            return None # Pass control back, note the finally clause is still executed 
        elif len(db_val) > 1:
            warnings.warn("More than 1 entry , warning :" 
                          "duplicate pseudo_prepids !")
        alignseqfileloc = db_val[-1][0] ## Get the last result
        vcf_loc = os.path.join(alignseqfileloc, "combined")
        vcf = os.path.join(vcf_loc, "{}.analysisReady.annotated.vcf".format(sample_name))
        if os.path.isfile(vcf):
            return vcf
        vcf += ".gz"
        if os.path.isfile(vcf):
            return vcf
        else:
            return None
    finally:
        if db.open:
            db.close()
   
class CreateGenomeBed(GATKFPipelineTask):
    """Use bedtools to create the input bed file for the coverage binning script"""
    def pre_shell_commands(self):
        self.shell_options.update(
            {"stdout":self.genome_cov_bed, "stderr":self.log_file})
        self.commands = [self.format_string(
            "{bedtools} genomecov -bga -ibam {recal_bam}")]

    def requires(self):
        return self.clone(PrintReads)

class CvgBinning(GATKFPipelineTask):
    def pre_shell_commands(self):
        try: # This syntax is needed to avoid race conditions in certain cases  
            if not os.path.isdir(self.cov_dir):
                os.makedirs(self.cov_dir)
        except OSError,e:
            if e.errno != 17:
                raise Exception("Problem Creating Directory : {0}".format(e))
            pass
        self.human_chromosomes = [str(x) for x in xrange(1, 23)] + ["X", "Y", "MT"]
        self.commands = [self.format_string(
            "{pypy} {coverage_binner} 1000 {name_prep} {genome_cov_bed} {cov_dir}")]

    def _run_post_success(self):
        os.remove(self.genome_cov_bed)

    def requires(self):
        return self.clone(CreateGenomeBed)

class GQBinning(GATKFPipelineTask):
    def pre_shell_commands(self):
        try:
            if not os.path.isdir(self.gq_dir):
                os.makedirs(self.gq_dir)
        except OSError,e:
            if e.errno != 17:
                raise Exception("Problem Creating Directory : {0}".format(e))
            pass
        self.human_chromosomes = [str(x) for x in xrange(1, 23)] + ["X", "Y", "MT"]
        self.binning_cmd = self.format_string(
            "{pypy} {gq_binner} 10000 {gvcf} {name_prep} {gq_dir}")
        self.commands = [self.binning_cmd]
        
    def requires(self):
        return self.clone(HaplotypeCaller)

class AlignmentMetrics(GATKFPipelineTask):
    def pre_shell_commands(self):
        self.alignment_metrics_raw = os.path.join(
            self.scratch_dir,
            "{sample_name}.{pseudo_prepid}.alignment.metrics.raw".format(
                sample_name=self.sample_name, pseudo_prepid=self.pseudo_prepid))
        cmd = self.format_string(
            "{java} -XX:ParallelGCThreads=4 -jar {picard} "
            "CollectAlignmentSummaryMetrics TMP_DIR={scratch_dir} "
            "VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE={ref_genome} "
            "INPUT={recal_bam} OUTPUT={alignment_metrics_raw} >& {log_file}")
        parser_cmd = self.format_string(
            "grep -v '^#' {alignment_metrics_raw} | awk -f {transpose_awk} "
            "> {alignment_metrics} 2> {log_file}")
        self.shell_options.update({"shell":True, "stdout":None, "stderr":None})
        self.commands = [cmd, parser_cmd]

    def requires(self):
        return self.clone(PrintReads)

class RunCvgMetrics(GATKFPipelineTask):
    def pre_shell_commands(self):
        """Run Picard CalculateHsMetrics,DepthOfCoverage or CollectWgsMetrics
        An optimization for the future is to split this up into multiple tasks since they are all independent
        of each other will speed things up significantly, IO can be an issue"""
        self.output_file = os.path.join(
            self.scratch_dir, self.name_prep + ".cvg.metrics.txt")
        self.raw_output_file = os.path.join(
            self.scratch_dir, self.sample_name + ".cvg.metrics.raw")
        self.output_file_ccds = os.path.join(
            self.scratch_dir, self.name_prep + ".cvg.metrics.ccds.txt")
        self.raw_output_file_ccds = os.path.join(
            self.scratch_dir, self.sample_name + ".cvg.metrics.ccds.raw")
        self.output_file_X = os.path.join(
            self.scratch_dir, self.name_prep + ".cvg.metrics.X.txt")
        self.raw_output_file_X = os.path.join(
            self.scratch_dir, self.sample_name + ".cvg.metrics.X.raw")
        self.output_file_Y = os.path.join(
            self.scratch_dir, self.name_prep + ".cvg.metrics.Y.txt")
        self.raw_output_file_Y = os.path.join(
            self.scratch_dir, self.sample_name + ".cvg.metrics.Y.raw")
        self.raw_output_file_cs = os.path.join(
            self.scratch_dir, self.sample_name + ".cvg.metrics.cs.raw")
        self.output_file_cs = os.path.join(
            self.scratch_dir, self.name_prep + ".cvg.metrics.cs.txt")
        self.DBID, self.prepID = getDBIDMaxPrepID(self.pseudo_prepid)
        if self.sample_type == "GENOME":
            ## Define shell commands to be run
            cvg_cmd = ("{java} -Xmx{max_mem}g -XX:ParallelGCThreads=4 -jar "
                            "{picard} CollectWgsMetrics VALIDATION_STRINGENCY=LENIENT "
                            "R={ref_genome} I={recal_bam} INTERVALS={file_target} "
                            "O={raw_output_file} MQ=20 Q=10 >> {log_file}")
            ## Run on the ccds regions only
            self.cvg_cmd1 = self.format_string(
                cvg_cmd, raw_output_file=self.raw_output_file_ccds,
                file_target=self.config_parameters["target_file"])
            ## Run across the genome
            self.cvg_cmd2 = self.format_string(
                "{java} -Xmx{max_mem}g -XX:ParallelGCThreads=4 -jar "
                "{picard} CollectWgsMetrics R={ref_genome} I={recal_bam} "
                "O={raw_output_file} MQ=20 Q=10 >> {log_file}")
            ## Run on X and Y Chromosomes only (across all regions there not just ccds)
            self.cvg_cmd3 = self.format_string(cvg_cmd,
                file_target=self.config_parameters["target_file_X"],
                raw_output_file=self.raw_output_file_X)
            self.cvg_cmd4 = self.format_string(cvg_cmd,
                file_target=self.config_parameters["target_file_Y"],
                raw_output_file=self.raw_output_file_Y)
        else:
            db = get_connection("seqdb")
            try:
                cur = db.cursor()
                query = ("SELECT region_file_lsrc FROM captureKit "
                         "INNER JOIN prepT ON prepT_name=prepT.exomeKit "
                         "WHERE prepID = {0} AND chr = 'X'".format(self.prepID))
                cur.execute(query)
                row = cur.fetchone()
                if row:
                    self.capture_file_X = row[0]
                elif self.sample_type == "CUSTOM_CAPTURE":
                    # we don't require X coverage for these
                    self.capture_file_X = None
                else:
                    raise ValueError("Couldn't find BED file for coverage on "
                                    "X chromosome")
                query = ("SELECT region_file_lsrc FROM captureKit "
                         "INNER JOIN prepT ON prepT_name=prepT.exomeKit "
                         "WHERE prepID = {0} AND chr = 'Y'".format(self.prepID))
                cur.execute(query)
                row = cur.fetchone()
                if row:
                    self.capture_file_Y = row[0]
                elif self.sample_type == "CUSTOM_CAPTURE":
                    self.capture_file_Y = None
                else:
                    raise ValueError("Couldn't find BED file for coverage on "
                                     "Y chromosome")
            finally:
                if db.open:
                    db.close()

            self.calc_capture_specificity = 1 ## We need to run an additional picard module for calculating capture specificity
            cvg_cmd = ("{java} -Xmx{max_mem}g -XX:ParallelGCThreads=4 -jar "
                       "{gatk} -T DepthOfCoverage -mbq 10 -mmq 20 "
                       "--omitIntervalStatistics --omitDepthOutputAtEachBase "
                       "--omitLocusTable -R {ref_genome} -I {recal_bam} "
                       "-ct 5 -ct 10 -ct 15 -ct 20 -L {file_target} "
                       "-o {raw_output_file} >> {log_file}")
            ## Run across ccds regions
            self.cvg_cmd1 = self.format_string(cvg_cmd,
                file_target=self.config_parameters["target_file"],
                raw_output_file=self.raw_output_file_ccds)
            ## Run across capture kit regions
            self.cvg_cmd2 = self.format_string(cvg_cmd,
                file_target=self.capture_kit_bed,
                raw_output_file=self.raw_output_file)
            ## Run on X and Y Chromosomes only (across all regions there not just ccds)
            if self.capture_file_X:
                self.cvg_cmd3 = self.format_string(cvg_cmd,
                    file_target=self.capture_file_X,
                    raw_output_file=self.raw_output_file_X)
            if self.capture_file_Y:
                self.cvg_cmd4 = self.format_string(cvg_cmd,
                    file_target=self.capture_file_Y,
                    raw_output_file=self.raw_output_file_Y)
            ## Run PicardHsMetrics for Capture Specificity
            self.output_bait_file = os.path.join(self.scratch_dir, "targets.interval")
            self.convert_cmd = self.format_string(
                "{java} -jar {picard} BedToIntervalList I={capture_kit_bed} "
                "SD={seqdict_file} OUTPUT={output_bait_file} >> {log_file}")
            self.cvg_cmd5 = self.format_string(
                "{java} -Xmx{max_mem}g -XX:ParallelGCThreads=4 -jar {picard} "
                "CollectHsMetrics BI={output_bait_file} TI={output_bait_file} "
                "VALIDATION_STRINGENCY=SILENT "
                "METRIC_ACCUMULATION_LEVEL=ALL_READS "
                "I={recal_bam} O={raw_output_file_cs} MQ=20 Q=10 >> {log_file}")
            ## DepthOfCoverage output has an extra suffix at the end
            self.raw_output_file_ccds = self.raw_output_file_ccds + '.sample_summary'
            self.raw_output_file = self.raw_output_file + '.sample_summary'
            self.raw_output_file_X = self.raw_output_file_X + '.sample_summary'
            self.raw_output_file_Y = self.raw_output_file_Y + '.sample_summary'

        ## An intermediate parsed file is created for extracting info for db update later, this remains the same for exomes and genomes
        parser_cmd = ("grep -v '^#' {raw_output_file} | "
                           "awk -f {transpose_awk} > {output_file}")
        self.parser_cmd1 = self.format_string(
            parser_cmd, raw_output_file=self.raw_output_file_ccds, output_file=self.output_file_ccds)
        self.parser_cmd2 = self.format_string(
            parser_cmd, raw_output_file=self.raw_output_file, output_file=self.output_file)
        if self.capture_file_X:
            self.parser_cmd3 = self.format_string(
                parser_cmd, raw_output_file=self.raw_output_file_X, output_file=self.output_file_X)
        if self.capture_file_Y:
            self.parser_cmd4 = self.format_string(
                parser_cmd, raw_output_file=self.raw_output_file_Y, output_file=self.output_file_Y)
        if self.sample_type != "GENOME": ## i.e. Exome or CustomCapture
            self.parser_cmd5 = self.format_string(
                parser_cmd, raw_output_file=self.raw_output_file_cs,
                output_file=self.output_file_cs)
        self.shell_options.update({"stdout":None, "stderr":None, "shell":True})
        self.commands = [
            self.cvg_cmd1, self.parser_cmd1, self.cvg_cmd2, self.parser_cmd2]
        if self.capture_file_X:
            self.commands.extend([self.cvg_cmd3, self.parser_cmd3])
        if self.capture_file_Y:
            self.commands.extend([self.cvg_cmd4, self.parser_cmd4])
        if self.sample_type != "GENOME":
            self.commands.extend(
                [self.convert_cmd, self.cvg_cmd5, self.parser_cmd5])

    def requires(self):
        yield self.clone(PrintReads)

class DuplicateMetrics(GATKFPipelineTask):
    """ Parse Duplicate Metrics from dragen logs """
    def pre_shell_commands(self):
        self.dragen_log = os.path.join(
            self.log_dir, self.name_prep + ".dragen.out")
        self.duplicates_file = os.path.join(
            self.scratch_dir, self.name_prep + ".duplicates.txt")
        if not os.path.isfile(self.dragen_log):
            raise Exception("The dragen log file could not be found!")
        perc_duplicates = None
        with open(self.dragen_log) as d:
            for line in d:
                line = line.strip()
                if line.endswith("%) duplicates marked"):
                    perc_duplicates = float(line.split("(")[-1].split("%")[0])
                    break
        if not perc_duplicates:
            # try alternate format
            with open(self.dragen_log) as d:
                for line in d:
                    line = line.strip()
                    if line.startswith("Number of duplicate reads"):
                        perc_duplicates = float(line.split("[")[-1].split("]")[0])
        if perc_duplicates:
            with open(self.duplicates_file, "w") as out:
                out.write(str(perc_duplicates) + "\n")
        else:
            raise ValueError("Could not find duplicate metrics in dragen log!")

class VariantCallingMetrics(GATKFPipelineTask):
    def pre_shell_commands(self):
        self.metrics_raw = os.path.join(
            self.scratch_dir, self.sample_name + ".raw")
        self.metrics_raw_summary = self.metrics_raw + ".variant_calling_summary_metrics"
        self.metrics_raw_detail = self.metrics_raw + "variant_calling_detail_metrics"
        self.metrics_file = os.path.join(
            self.scratch_dir, self.name_prep + ".variant_calling_summary_metrics.txt")
        self.metrics_cmd = self.format_string(
            "{java} -XX:ParallelGCThreads=4 -jar {picard} CollectVariantCallingMetrics "
            "INPUT={annotated_vcf_gz} OUTPUT={metrics_raw} DBSNP={dbSNP} >& {log_file}")
        self.parser_cmd = self.format_string(
            "grep -v '^#' {metrics_raw_summary} | awk -f {transpose_awk} "
            "> {metrics_file}")
        self.shell_options.update(
            {"stdout":None, "stderr":None, "shell":True})
        self.commands = [self.metrics_cmd, self.parser_cmd]

    def requires(self):
        return self.clone(AnnotateVCF)

class GenotypeConcordance(GATKFPipelineTask):
    """ Run the gatk tool for evaluation of the variant calls
    with older production pipeline variant calls """
    n_cpu = 7
    def pre_shell_commands(self):
        self.truth_vcf = os.path.join(self.scratch_dir, "truth.vcf")
        self.eval_vcf = os.path.join(self.scratch_dir, "eval.vcf")
        self.concordance_metrics = os.path.join(
            self.scratch_dir, self.name_prep + ".genotype_concordance_metrics.txt")
        self.concordance_cmd = self.format_string(
            "{java} -XX:ParallelGCThreads=4 -jar {gatk} -T GenotypeConcordance "
            "-R {ref_genome} -eval {eval_vcf} -comp {truth_vcf} "
            "-o {concordance_metrics} -U ALLOW_SEQ_DICT_INCOMPATIBILITY")
        production_vcf = get_productionvcf(
            self.pseudo_prepid, self.sample_name,self.sample_type)
        # wait for Avere to cooperate before reading the production VCF
        while True:
            if check_avere():
                break
            else:
                sleep(25 + 10 * random())
        self.subset_vcf(production_vcf, self.truth_vcf)
        self.subset_vcf(self.annotated_vcf_gz, self.eval_vcf)
        for idx in [vcf +  ".idx" for vcf in (self.truth_vcf, self.eval_vcf)]:
            if os.path.isfile(idx):
                os.remove(idx)
        self.commands = [self.format_string(self.concordance_cmd)]

    def subset_vcf(self, vcf_fn, vcf_out_fn):
        """ Get only PASS variants based on the FILTER field
        Returns : Path to the output vcf
        """
        with get_fh(vcf_fn) as vcf_in, open(vcf_out_fn, "w") as vcf_out:
            for line in vcf_in:
                vcf_out.write(line)
                if line.startswith("#CHROM"):
                    FILTER_idx = fields.split("\t").index("FILTER")
                    break
            for line in vcf_in:
                if line.split("\t")[FILTER_idx] == "PASS":
                    vcf_out.write(line)

    def requires(self):
        return self.clone(CombineVariants)

class ContaminationCheck(GATKFPipelineTask):
    """ Run VerifyBamID to check for sample contamination """
    def pre_shell_commands(self):
        self.contamination_raw = os.path.join(
            self.scratch_dir, self.sample_name + ".contamination.raw")
        self.contamination_final = os.path.join(
            self.scratch_dir, self.name_prep + ".contamination.selfSM.txt")
        self.contamination_cmd = self.format_string(
            "{verifybamid} --vcf {contam1000g_vcf} --bam {recal_bam} --out "
            "{contamination_raw} --verbose --ignoreRG --maxDepth 1000 --precise "
            "> {log_file} 2> {err}")
        self.parser_cmd = self.format_string(
            "awk -f {transpose_awk} {contamination_raw}.selfSM > {contamination_final}")
        self.shell_options.update(
            {"stdout":None, "stderr":None, "shell":True})
        self.commands = [self.contamination_cmd, self.parser_cmd]

    def requires(self):
        return self.clone(PrintReads)

def update_alignseqfile(cur, location, pseudo_prepid):
    update_statement = """
    UPDATE dragen_qc_metrics SET AlignSeqFileLoc = "{location}"
    WHERE pseudo_prepid = {pseudo_prepid}""".format(
        location=location, pseudo_prepid=pseudo_prepid)
    try:
        cur.execute(update_statement)
    except MySQLdb.Error, e:
        raise Exception("ERROR %d IN CONNECTION: %s" % (e.args[0], e.args[1]))

class UpdateSeqdbMetrics(GATKFPipelineTask):
    def pre_shell_commands(self):
        ## Generic query to be used for updates 
        self.update_statement = """
        UPDATE {table} SET {field} = "{value}"
        WHERE CHGVID = "{sample_name}" AND seqType = "{sample_type}"
            AND pseudo_prepid = "{pseudo_prepid}" """
        self.query_statement = """
        SELECT {field}
        FROM {table}
        WHERE CHGVID = '{sample_name}' AND seqType = '{sample_type}'
            AND prepid = {prep_id}"""
        self.issue_contamination_warning = False
        ## The output files from the tasks run 
        self.cvg_out = os.path.join(
            self.scratch_dir, self.name_prep + ".cvg.metrics.txt")
        self.cvg_ccds_out = os.path.join(
            self.scratch_dir, self.name_prep + ".cvg.metrics.ccds.txt")
        self.cvg_X_out = os.path.join(
            self.scratch_dir, self.name_prep + ".cvg.metrics.X.txt")
        self.cvg_Y_out = os.path.join(
            self.scratch_dir, self.name_prep + ".cvg.metrics.Y.txt")
        self.cvg_cs_out = os.path.join(
            self.scratch_dir, self.name_prep + ".cvg.metrics.cs.txt")
        self.dup_out = os.path.join(
            self.scratch_dir, self.name_prep + ".duplicates.txt")
        self.variant_call_out = os.path.join(
            self.scratch_dir, self.name_prep + ".variant_calling_summary_metrics.txt")
        self.geno_concordance_out = os.path.join(
            self.scratch_dir, self.name_prep + ".genotype_concordance_metrics.txt")
        self.contamination_out = os.path.join(
            self.scratch_dir, self.name_prep + ".contamination.selfSM.txt")
        ## The qc table to update
        self.qc_table = "dragen_qc_metrics"
        ## The new lean equivalent of seqdbClone 
        self.master_table = "seqdbClone"
        self.DBID,self.prepID = getDBIDMaxPrepID(self.pseudo_prepid)
        raw_vcfs = os.path.join(
            self.scratch_dir, self.name_prep + ".raw*vcf*")
        indel_vcfs = os.path.join(
            self.scratch_dir, self.name_prep + ".indel*")
        snp_vcfs = os.path.join(
            self.scratch_dir, self.name_prep + ".snp*")
        analysis_vcfs = os.path.join(
            self.scratch_dir, self.name_prep + ".analysisReady.vcf*")
        fixed_vcf = os.path.join(
            self.scratch_dir, self.name_prep + ".analysisReady.fixed.vcf")
        tmp_vcf = os.path.join(
            self.scratch_dir, self.name_prep + ".tmp.vcf")
        self.files_to_remove = [raw_vcfs, indel_vcfs, snp_vcfs, analysis_vcfs, fixed_vcf, tmp_vcf]
        self.cvg_parse = {"EXOME":{
            "all":{"mean":qc_metrics().mean_cvg, "granular_median":qc_metrics().median_cvg,
                   "%_bases_above_5":qc_metrics().pct_bases5X, "%_bases_above_10":qc_metrics().pct_bases10X,
                   "%_bases_above_15":qc_metrics().pct_bases15X, "%_bases_above_20":qc_metrics().pct_bases20X},
            "X":{"mean":qc_metrics().mean_X_cvg},
            "Y":{"mean":qc_metrics().mean_Y_cvg},
            "ccds":{"mean":qc_metrics().mean_ccds_cvg, "%_bases_above_5":qc_metrics().pct_ccds_bases5X,
                    "%_bases_above_10":qc_metrics().pct_ccds_bases10X,
                    "%_bases_above_15":qc_metrics().pct_ccds_bases15X,
                    "%_bases_above_20":qc_metrics().pct_ccds_bases20X},
            "cs":{"ON_BAIT_VS_SELECTED":qc_metrics().capture_specificity}}, "GENOME":{
            "all":{"MEAN_COVERAGE":qc_metrics().mean_cvg, "MEDIAN_COVERAGE":qc_metrics().median_cvg,
                   "PCT_5X":qc_metrics().pct_bases5X, "PCT_10X":qc_metrics().pct_bases10X,
                   "PCT_15X":qc_metrics().pct_bases15X, "PCT_20X":qc_metrics().pct_bases20X},
            "X":{"MEAN_COVERAGE":qc_metrics().mean_X_cvg},
            "Y":{"MEAN_COVERAGE":qc_metrics().mean_Y_cvg},
            "ccds":{"MEAN_COVERAGE":qc_metrics().mean_cvg, "PCT_5X":qc_metrics().pct_bases5X,
                    "PCT_10X":qc_metrics().pct_bases10X, "PCT_15X":qc_metrics().pct_bases15X,
                    "PCT_20X":qc_metrics().pct_bases20X}}}
        try:
            self.db = get_connection("seqdb")
            self.cur = self.db.cursor()
            self.get_variant_counts()
            self.add_sample_to_qc()
            self.update_alignment_metrics()
            self.update_coverage_metrics('all')
            self.update_coverage_metrics('ccds')
            self.update_coverage_metrics('X')
            self.update_coverage_metrics('Y')
            if self.sample_type != 'GENOME':
                self.update_coverage_metrics('cs')
            self.update_duplicates()
            self.update_variant_calling_metrics()
            if get_productionvcf(
                self.pseudo_prepid, self.sample_name, self.sample_type) is not None:
                self.update_genotype_concordance_metrics()
            self.update_contamination_metrics()
            self.update_seqgender()
            self.update_qc_message()
            ## Update alignseqfile location here to keep track of scratch location
            # Due to bad decision the data location does not include the
            # sample's name and prep ID
            update_alignseqfile(
                self.cur, os.path.dirname(self.scratch_dir), self.pseudo_prepid)
            self.db.commit()
            for tmp_file in self.files_to_remove:
                for fn in glob(tmp_file):
                    os.remove(fn)
        finally:
            if self.db.open:
                self.db.close()

    def add_sample_to_qc(self):
        """
        Initialize the sample in the qc table, use update statements later to update qc metrics
        """
        query = self.format_string("""
        INSERT INTO {qc_table} (CHGVID, seqType, pseudo_prepid)
        VALUES ("{sample_name}", "{sample_type}", {pseudo_prepid})""")
        self.cur.execute(query)

    def check_qc(self):
        """
        Return any failed QC checks and associated messages
        """
        qc_failures = OrderedDict()
        for qc, msg in (("check_alignment_rate", "Failed Alignment Check"),
                        ("check_duplicates", "Failed Duplicate Check"),
                        ("check_variant_calling", "Failed VariantCallingCheck"),
                        ("check_coverage", "Failed Coverage Check"),
                        ("check_contamination", "Failed Contamination Check")):
            if not getattr(self, qc)():
                qc_failures[qc] = msg
        return qc_failures

    def update_qc_message(self):
        """ Update qc message in statusT based on qc check, will call the pass_qc
        function to perform the check """
        self.failed_qc = "QC review needed"
        self.passed_qc = "Passed Bioinfo QC"
        message = []
        self.db = get_connection("seqdb")
        qc_failures = self.check_qc()
        if qc_failures:
            message.append(self.failed_qc)
            if self.issue_contamination_warning == True: ## Check for contamination warning
                message.append("Warning sample contamination is high, but below qc fail threshold")
            message.extend(qc_failures.values())
        else:
            message.append(self.passed_qc)
        final_message = ';'.join(message)
        self.update_database(self.qc_table,'QCMessage',final_message)

    def update_database(self, table, field, value):
        self.cur.execute(self.format_string(
            self.update_statement, table=table, field=field, value=value))

    def get_metrics(self, query):
        try:
            self.cur.execute(query)
            db_val = self.cur.fetchall()
            return db_val
        except MySQLdb.Error, e:
            raise ValueError(query)

    def update_alignment_metrics(self):
        self.alignment_metrics_map = {"TOTAL_READS":qc_metrics().total_reads,
                                      "PCT_PF_READS_ALIGNED":qc_metrics().pct_reads_aligned,
                                      "PF_MISMATCH_RATE":qc_metrics().pct_mismatch_rate}

        with open(self.alignment_metrics) as alignment_metrics_fh:
            for line in alignment_metrics_fh:
                contents = line.strip().split(" ")
                if len(contents) == 4:
                    field = contents[0]
                    value = contents[-1]
                    if field in self.alignment_metrics_map:
                        db_field = self.alignment_metrics_map[field]
                        self.update_database(self.qc_table, db_field, value)

    def update_coverage_metrics(self,file_type):
        if file_type == "ccds":
            metrics_file = self.cvg_ccds_out
        elif file_type == "X":
            if (self.sample_type == "CUSTOM_CAPTURE"
                and not os.path.isfile(self.cvg_X_out)):
                return
            else:
                metrics_file = self.cvg_X_out
        elif file_type == "Y":
            if (self.sample_type == "CUSTOM_CAPTURE"
                and not os.path.isfile(self.cvg_Y_out)):
                return
            else:
                metrics_file = self.cvg_Y_out
        elif file_type == "cs":
            metrics_file = self.cvg_cs_out
        elif file_type == "all":
            metrics_file = self.cvg_out
        sample_type = "GENOME" if self.sample_type == "GENOME" else "EXOME"
        metrics_hash = self.cvg_parse[sample_type][file_type]

        with open(metrics_file) as metrics_fh:
            for line in metrics_fh:
                contents = line.strip().split(" ")
                if len(contents) > 1:
                    field, value  = contents[:2]
                    if field in metrics_hash:
                        db_field = metrics_hash[field]
                        self.update_database(self.qc_table, db_field, value)
                        if file_type == "all":
                            if db_field == "OverallCoverage":
                                mean_cvg = float(value)
                            if db_field == "MedianCoverage":
                                median_cvg = float(value)

        if file_type == "all":
            db_field = qc_metrics().mean_median_cvg
            value = mean_cvg / median_cvg
            self.update_database(self.qc_table, db_field, value)

    def update_duplicates(self):
        with open(self.dup_out) as duplicates_fh:
            pct_duplicates = float(duplicates_fh.next())
        self.update_database(
            self.qc_table, qc_metrics().pct_duplicate_reads, pct_duplicates)

    def update_variant_calling_metrics(self):
        ## Need to calculate overall titv since that is not available in the output file
        variant_call_parse = {
            "TOTAL_SNPS":qc_metrics().total_snps, "PCT_DBSNP":qc_metrics().pct_dbsnp_snps,
            "TOTAL_INDELS":qc_metrics().total_indels, "PCT_DBSNP_INDELS":qc_metrics().pct_dbsnp_indels}
        temp = {}
        flag = 0 
        with open(self.variant_call_out) as variant_call_metrics_fh:
            for line in variant_call_metrics_fh:
                field, value = line.strip().split(" ")[:2]
                if field in variant_call_parse:
                    db_field = variant_call_parse[field]
                    self.update_database(self.qc_table, db_field, value)
                temp[field] = float(value)

        overall_titv = (
            (temp["NOVEL_SNPS"] * temp["NOVEL_TITV"] + temp["NUM_IN_DB_SNP"] * temp["DBSNP_TITV"]) /
            temp["TOTAL_SNPS"])
        self.update_database(self.qc_table, qc_metrics().titv, overall_titv)

        homhet_ratio = (float(self.all_snv_hom + self.all_indel_hom) /
                        (self.all_snv_het + self.all_indel_het))
        self.update_database(self.qc_table, qc_metrics().homhet_ratio, homhet_ratio)

        snv_homhet_ratio = float(self.all_snv_hom) / self.all_snv_het
        self.update_database(self.qc_table, qc_metrics().snv_homhet_ratio, snv_homhet_ratio)

        indel_homhet_ratio = float(self.all_indel_hom) / self.all_indel_het
        self.update_database(self.qc_table, qc_metrics().indel_homhet_ratio, indel_homhet_ratio)

        x_homhet_ratio = (float(self.X_snv_hom + self.X_indel_hom) /
                          (self.X_snv_het + self.X_indel_het))
        self.update_database(self.qc_table, qc_metrics().x_homhet_ratio, x_homhet_ratio)
             
    def update_genotype_concordance_metrics(self):
        genotype_concordance = {
            "ALLELES_MATCH":1, "ALLELES_DO_NOT_MATCH":1, "EVAL_ONLY":1, "TRUTH_ONLY":1}
        self.temp_geno_concordance = os.path.join(
            self.scratch_dir, "temp_geno_concordance.txt")
        cmd = self.format_string(
            "grep -A 2 '#:GATKTable:SiteConcordance_Summary:Site-level summary statistics' "
            "{geno_concordance_out} | grep -v '#' | awk -f {transpose_awk} > {temp_geno_concordance}")
        proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        proc.wait()
        if proc.returncode:
            print >> self.LOG_FILE,subprocess.CalledProcessError(proc.returncode,self.get_het_cmd)
        with open(self.temp_geno_concordance,'r') as concordance_fh:
            for line in concordance_fh:
                field, value = line.strip().split(' ')[:2]
                if field in genotype_concordance:
                    genotype_concordance[field] = int(value)

        concordance = (
            float(genotype_concordance['ALLELES_MATCH']) /
            (genotype_concordance['ALLELES_DO_NOT_MATCH'] + genotype_concordance['ALLELES_MATCH']
             + genotype_concordance['EVAL_ONLY'] + genotype_concordance['TRUTH_ONLY']))
        self.update_database(self.qc_table, qc_metrics().concordance, concordance)

    def update_contamination_metrics(self):
        with open(self.contamination_out) as contamination_fh:
            for line in contamination_fh:
                field, value = line.strip().split(' ')[:2]
                if field == "FREEMIX":
                    db_field = qc_metrics().contamination_value
                    self.update_database(self.qc_table, db_field,value)
                    return

    def update_seqgender(self):
        """
        Calculates the X/Y coverage for Exomes and Genomes and updates the seq gender
        based on existing rules : https://redmine.igm.cumc.columbia.edu/projects/biopipeline/wiki/Gender_checks
        For custom capture regions the rules are more complicated. 
        Also updates the database with het/hom ratio on X chromosome, the num_homX and the num_hetX
        Args : Nothing
        Returns : String; One of the following values : [Male,Female,Ambiguous]
        """

        if self.sample_type in ("EXOME","GENOME"):
            query = """
            SELECT {cvg_type}
            FROM {qc_table}
            WHERE CHGVID = "{sample_name}" AND seqType = "{sample_type}"
                AND pseudo_prepid = {pseudo_prepid}"""
            result = self.get_metrics(
                self.format_string(query, cvg_type=qc_metrics().mean_X_cvg))
            if result:
                X_cvg = float(result[0][0])
                result = self.get_metrics(
                    self.format_string(query, cvg_type=qc_metrics().mean_Y_cvg))
                if result:
                    Y_cvg = float(result[0][0])
                    if not isclose(0.0, Y_cvg):
                        ratio = X_cvg / Y_cvg
                        self.update_database(
                            self.qc_table, qc_metrics().xy_avg_ratio, ratio)
                    if isclose(0.0, Y_cvg) or ratio > 5:
                        seq_gender = "F"
                    elif ratio < 2:
                        seq_gender = "M"
                    else:
                        seq_gender = "Ambiguous"
                else:
                    seq_gender = "Ambiguous"
            else:
                seq_gender = "Ambiguous"

        else : ## Custom Capture Samples
            het = self.X_snv_het + self.X_indel_het
            hom = self.X_snv_hom + self.X_indel_hom
            if hom:
                ratio = float(het) / hom
            if not (het or hom) or (0.26 <= ratio <= 0.58):
                seq_gender = "Ambiguous"
            elif not het:
                seq_gender = "M"
            elif not hom:
                seq_gender = "F"
            elif ratio < 0.26:
                seq_gender = "M"
            elif ratio > 0.58:
                seq_gender = "F"

        self.update_database(self.qc_table, 'SeqGender', seq_gender)

    def check_alignment_rate(self):
        """
        Checks if the alignment metrics matches the appropriate thresholds
        Currently using only the perc_reads_aligned
        """

        query = self.format_string(
            """SELECT {pct_reads_aligned}
            FROM {qc_table}
            WHERE CHGVID = "{sample_name}" AND seqType = "{sample_type}" AND
                pseudo_prepid = {pseudo_prepid}""")
        result = self.get_metrics(query)
        perc_reads_aligned = float(result[0][0])
        return (perc_reads_aligned >= 0.70)

    def check_coverage(self):
        """
        Checks if the coverage metrics matches the appropriate thresholds
        """

        query = self.format_string(
            """SELECT {pct_bases5X}, {mean_cvg}
            FROM {qc_table}
            WHERE CHGVID = "{sample_name}" AND seqType = "{sample_type}" AND
                pseudo_prepid = {pseudo_prepid}""")
        result = self.get_metrics(query)
        pct_bases_5X = float(result[0][0])
        mean_cvg = float(result[0][1])
        return (pct_bases_5X >= 0.9 and mean_cvg >= 25)

    def check_duplicates(self):
        """
        Checks if the duplicate metrics matches the appropriate thresholds
        Currently using only the perc_duplicate_reads
        """
        query = self.format_string(
            """SELECT {pct_duplicate_reads}
            FROM {qc_table}
            WHERE CHGVID = "{sample_name}" AND seqType = "{sample_type}" AND
                pseudo_prepid = {pseudo_prepid}""")
        result = self.get_metrics(query)
        perc_duplicate_reads = float(result[0][0])
        if self.sample_type == "GENOME":
            return (perc_duplicate_reads <= 20)

        else:
            return (perc_duplicate_reads <= 30)

    def check_variant_calling(self):
        """
        Check if a minimal number of SNVs have been called, per sequencing type
        Could very well fail on custom capture
        """
        query = self.format_string(
            """SELECT {total_snps}
            FROM {qc_table} WHERE CHGVID = "{sample_name}" AND seqType = "{sample_type}"
            AND pseudo_prepid = {pseudo_prepid}""")
        result = self.get_metrics(query)
        num_snvs = int(result[0][0])

        if self.sample_type == "GENOME":
            return num_snvs > 3000000
        else:
            return num_snvs > 100000
        # arguably should use something else for CUSTOM_CAPTURE...

    def check_contamination(self):
        """
        Checks if the contamination metrics matches the appropriate thresholds
        """
        query = self.format_string(
            """SELECT {contamination_value}
            FROM {qc_table}
            WHERE CHGVID = "{sample_name}" AND seqType = "{sample_type}"
                AND pseudo_prepid = {pseudo_prepid}""")
        result = self.get_metrics(query)
        contamination_value = float(result[0][0])
        if contamination_value < 0.08:
            if contamination_value >= 0.05:
                self.issue_contamination_warning = True
            return True
        elif contamination_value >= 0.08:
            return False

    def check_gender(self):
        """
        Query seqdbClone and check whether the seq and selfdecl gender matchup
        """
        base_query = """SELECT {sex_field}
        FROM {master_table}
        WHERE CHGVID = "{sample_name}" AND seqType = "{sample_type}"
        AND pseudo_prepid = {pseudo_prepid}"""
        declared_query = self.format_string(base_query, sex_field="SelfDeclGender")
        result = self.get_metrics(declared_query)[0][0]
        sequenced_query = self.format_string(base_query, sex_field="SeqGender")
        seq_sex = self.get_metrics(sequenced_query)[0][0]
        return self_declared_sex == seq_sex

    def get_variant_counts(self):
        """Set the number of het & hom variant counts, divided up by X vs.
        others
        Additionally require X to have MQ >= 40, QD > 1, and GQ >= 20 to
        contribute to X count
        """
        variants = {"X":{"snv":Counter(), "indel":Counter()},
                    "all":{"snv":Counter(), "indel":Counter()}}
        with gopen(self.annotated_vcf_gz) as vcf_fh:
            for line in vcf_fh:
                if line.startswith("#CHROM"):
                    break
            header_fields = line[1:].strip().split("\t")
            for line in vcf_fh:
                fields = dict(zip(header_fields, line.strip().split("\t")))
                lREF = len(fields["REF"])
                multiallelic = "," in fields["ALT"]
                if multiallelic:
                    nsnvs = len([1 for alt in fields["ALT"].split(",")[:2]
                                 if lREF == len(alt)])
                    # they're het at a multiallelic site
                    variants["all"]["snv"]["het"] += nsnvs
                    variants["all"]["indel"]["het"] += (2 - nsnvs)
                else:
                    variant_type = (
                        "snv" if lREF == len(fields["ALT"]) else "indel")
                    call_dict = dict(
                        zip(fields["FORMAT"].split(":"),
                            fields[self.sample_name].split(":")))
                    if len(set(call_dict["GT"].split("/"))) == 1:
                        call = "hom"
                    else:
                        call = "het"
                    variants["all"][variant_type][call] += 1
                if fields["CHROM"] == "X":
                    INFO_dict = dict([field.split("=", 1) for field in
                                      fields["INFO"].split(";") if "=" in field])
                    if (("MQ" in INFO_dict and int(INFO_dict["MQ"].split(".")[0]) > 39) and
                        ("QD" in INFO_dict and int(INFO_dict["QD"].split(".")[0]) > 1)):
                        if "GQ" in call_dict and int(call_dict["GQ"]) > 19:
                            if multiallelic:
                                variants["X"]["snv"]["het"] += nsnvs
                                variants["X"]["indel"]["het"] += (2 - nsnvs)
                            else:
                                variants["X"][variant_type][call] += 1
            self.X_snv_het = variants["X"]["snv"]["het"]
            self.X_snv_hom = variants["X"]["snv"]["hom"]
            self.X_indel_het = variants["X"]["indel"]["het"]
            self.X_indel_hom = variants["X"]["indel"]["hom"]
            self.all_snv_het = variants["all"]["snv"]["het"]
            self.all_snv_hom = variants["all"]["snv"]["hom"]
            self.all_indel_het = variants["all"]["indel"]["het"]
            self.all_indel_hom = variants["all"]["indel"]["hom"]

    def requires(self):
        yield self.clone(GQBinning)
        yield self.clone(CvgBinning)
        yield self.clone(AlignmentMetrics)
        yield self.clone(RunCvgMetrics)
        yield self.clone(DuplicateMetrics)
        yield self.clone(VariantCallingMetrics)
        yield self.clone(ContaminationCheck)
        if get_productionvcf(
            self.pseudo_prepid,self.sample_name,self.sample_type) is not None:
            yield self.clone(GenotypeConcordance)
