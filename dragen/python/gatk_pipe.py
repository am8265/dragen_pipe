#!/nfs/goldstein/software/python2.7.7/bin/python2.7

import os
# import math
import sys
import subprocess
import luigi
from luigi.contrib.sge import SGEJobTask
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
from collections import Counter, defaultdict
from getpass import getuser
from pwd import getpwuid
from gzip import open as gopen
from string import Formatter
from math import ceil
from check_vcf import check_vcf
from check_bam import check_bam
from dragen_db_statements import *
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.realpath(__file__)))), "import", "src"))
from waldb_globals import *
from inspect import currentframe, getframeinfo
import re
import datetime

warnings.filterwarnings("error", category=MySQLdb.Warning)

class ArchiveDirectoryAlreadyExists(Exception):
    pass

class VCFCheckException(Exception):
    pass

class BAMCheckException(Exception):
    pass

class RsyncException(Exception):
    pass

class QCError(Exception):
    pass

def owner(fn):
    """Return the user name of the owner of the specified file
    """
    if os.path.exists(fn):
        return getpwuid(os.stat(fn).st_uid).pw_name
    else:
        raise OSError("{fn} does not exist!".format(fn=fn))

class programs(luigi.Config):
    """ Configuration class for instantiating parameters for this pipeline
    the values are read from luigi.cfg in the current folder """
    bgzip = luigi.InputFileParameter()
    gatk = luigi.InputFileParameter()
    java = luigi.InputFileParameter()
    picard = luigi.InputFileParameter()
    pypy = luigi.InputFileParameter()
    tabix = luigi.InputFileParameter()
    verifybamid = luigi.InputFileParameter()
    clineff = luigi.InputFileParameter()

class locations(luigi.Config):
    base = luigi.Parameter()

class pipeline_files(luigi.Config):
    gq_binner = luigi.InputFileParameter()
    clineff_cfg = luigi.InputFileParameter()
    target_file = luigi.InputFileParameter()
    target_file_X = luigi.InputFileParameter()
    target_file_Y = luigi.InputFileParameter()
    transpose_awk = luigi.InputFileParameter()
    def __init__(self, *args, **kwargs):
        if "DEBUG_INTERVALS" in os.environ:
            kwargs["target_file"] = os.getenv("DEBUG_INTERVALS")
        super(pipeline_files, self).__init__(*args, **kwargs)

class gatk_resources(luigi.Config):
    contam1000g_vcf = luigi.InputFileParameter()
    ref_genome = luigi.InputFileParameter()
    seqdict_file = luigi.InputFileParameter()
    interval = luigi.InputFileParameter()
    silly_arg = luigi.Parameter(default="")
    g1000 = luigi.InputFileParameter()
    hapmap = luigi.InputFileParameter()
    Mills1000g = luigi.InputFileParameter()
    omni = luigi.InputFileParameter()
    dbSNP = luigi.InputFileParameter()
    annotatedbSNP = luigi.InputFileParameter()
    exac = luigi.InputFileParameter()
    roche_bed = luigi.InputFileParameter()

    def __init__(self, *args, **kwargs):
        if "DEBUG_INTERVALS" in os.environ:
            kwargs["interval"] = os.getenv("DEBUG_INTERVALS")
            kwargs["silly_arg"] = " -L " + kwargs["interval"]
            if "DEBUG_DBSNP" in os.environ:
                kwargs["dbSNP"] = os.getenv("DEBUG_DBSNP")
        super(gatk_resources, self).__init__(*args, **kwargs)

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

class GATKFPipelineTask(GATKPipelineTask):
    task_name_format = "{task_family}.{sample_name}.{pseudo_prepid}"
    poll_time = luigi.IntParameter(
        default=120,
        description="the number of seconds to wait before rerunning qstat for a task")
    subdirectory = luigi.Parameter(
        significant=False, default="{sample_name}.{pseudo_prepid}",
        description="output to this (relative) path instead of a randomly generated directory name")
    prept_start_message = "GATK Pipeline" # any inherited class will update
    # prepT's status column with having started or failed the GATK pipeline as
    # appropriate
    dont_remove_tmp_dir = True
    # don't delete temp dirs ever for any Task other than ArchiveSample,
    # as the directories are shared across Tasks

    def arseholes(self):

        if self.pipeline_step_id > 12: # if self.pipeline_step_id >= 23:

            check='/nfs/goldstein/software/tabix-0.2.6/tabix {} 22 | head -10 | wc -l'.format(self.gvcf)

            chr22check = subprocess.Popen(check,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            count = chr22check.stdout.readline().rstrip();

            if count == '10': # f'ing lazy
            # if chr22check.stdout.readline().rstrip() == '10': # f'ing lazy
                print ('g.vcf seems fine')
            else:
                db = get_connection("seqdb")
                cur = db.cursor()
                wtf = "update prepT set status = 'gVCF Truncated', status_time = UNIX_TIMESTAMP(CURRENT_TIMESTAMP()) where experiment_id = {}".format(self.pseudo_prepid)
                cur.execute(wtf)
                wtf = "update dragen_sample_metadata set is_merged = 90124 where experiment_id = {}".format(self.pseudo_prepid)
                cur.execute(wtf)
                db.commit()
                os._exit(1) 

        if self.pipeline_step_id > 23:

            check='/nfs/goldstein/software/tabix-0.2.6/tabix {}/{}.{}.analysisReady.vcf.gz 22 | head -10 | wc -l'.format(
              self.scratch_dir,self.sample_name,self.pseudo_prepid
            )
            chr22check = subprocess.Popen(check,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            if chr22check.stdout.readline().rstrip()=="10":
                print("analysis vcf seems fine");
            else:
                db = get_connection("seqdb")
                cur = db.cursor()
                wtf = "update prepT set status = 'VCF Truncated', status_time = UNIX_TIMESTAMP(CURRENT_TIMESTAMP()) where experiment_id = {}".format(self.pseudo_prepid)
                cur.execute(wtf)
                wtf = "update dragen_sample_metadata set is_merged = 90234 where experiment_id = {}".format(self.pseudo_prepid)
                cur.execute(wtf)
                db.commit()
                os._exit(1) 

        if self.pipeline_step_id < 14:
            return

        metrics = os.path.join(self.scratch_dir, self.name_prep + ".cvg.metrics.txt")
        mean=-0.0

        with open(metrics,"r") as f:
            for l in f:
                s=l.split(' ') 
                if s[0] == 'mean':
                    meany=float(s[1])
                    if meany==float(s[2]):
                        mean=meany

                    break
                else:
                    continue 
        if mean < 23.0:

            db = get_connection("seqdb")
            cur = db.cursor()

            # try:
            wtf = "update dragen_sample_metadata set is_merged = -2 where experiment_id = {}".format(self.pseudo_prepid)
            cur.execute(wtf)
            wtf = "update prepT set is_released = 0, status = 'Failed/Low-Qual Sample; Should never have been released ({}x) - will require deprecation', status_time = UNIX_TIMESTAMP(CURRENT_TIMESTAMP()) where experiment_id = {}".format(mean,self.pseudo_prepid)
            cur.execute(wtf)
            wtf = "update Experiment set is_released = 'release_rejected' where id = {}".format(self.pseudo_prepid)
            try:
                cur.execute(wtf)
            except:
                print("this is prolly cos there's nothing to change")

            db.commit()

            os._exit(1)

            # raise Exception("\n\nthis sample shouldn't be in the pipeline at all - mean cov = {}\n\n".format(mean))
            # exit(1)

            # cur.execute("update prepT set status = '{}', status_time = CURRENT_TIMESTAMP() where p_prepid = {}".format(status,self.pseudo_prepid))
            # except MySQLdb.Error, e:
                # raise Exception("ERROR %d IN CONNECTION: %s" % (e.args[0], e.args[1]))
            # finally:
                # if db.open:
                    # db.close()

        else:
            print("YAY")

    def __init__(self, *args, **kwargs):
        self.config_parameters = {}
        for cls in (programs, locations, pipeline_files, gatk_resources, qc_metrics):
            self.config_parameters.update(cls().__dict__)
        if "DEBUG_INTERVALS" in os.environ:
            kwargs["poll_time"] = 10 # use a shorter qstat poll time if running debug mode
            if "DEBUG_BED" in os.environ:
                kwargs["capture_kit_bed"] = os.getenv("DEBUG_BED")
        if self.sample_type == "GENOME_AS_FAKE_EXOME":
            self.config_parameters["interval"] = self.config_parameters["roche_bed"]
        super(GATKFPipelineTask, self).__init__(*args, **kwargs)
        self.config_parameters.update(self.__dict__)

    def format_string(self, s, **kwargs):
        """Format the string s with any parameters in config_parameters, self, and kwargs
        (priority if overlap is given to kwargs then self)
        """
        for field_name in (record[1] for record in list(Formatter().parse(s))):
            if (field_name and field_name not in self.config_parameters
                and hasattr(self, field_name)):
                # not all attributes get set in __init__ such that they appear
                # in __dict__ so this is a workaround to get them
                self.config_parameters[field_name] = getattr(self, field_name)
        self.config_parameters.update(kwargs)
        return s.format(**self.config_parameters)

    def archive_helper(self):

        self.data_to_copy = [ "{pipeline_tarball}", "{recal_bam}", "{recal_bam_index}", "{annotated_vcf_gz}", 
          "{annotated_vcf_gz_index}", "{gvcf}", "{gvcf_index}", "{cvg_tarball}", "{gq_tarball}", "{raw_coverage}"
        ]

        if self.sample_type == "GENOME":
            if os.path.isfile(self.original_vcf_gz):
                self.data_to_copy.append("{original_vcf_gz}")

        # copy original BAM in this case as well in case we want to reprocess as an actual genome laster
        if self.sample_type == "GENOME_AS_FAKE_EXOME":
            self.data_to_copy.extend(["{scratch_bam}", "{scratch_bam}.bai"])

    def update_sample_status(self, status):

        status = 'Pipeline ('+status+')'
        db = get_connection("seqdb")
        try:
            cur = db.cursor()
            try:
                wtf = "update prepT set status = '{}', status_time = UNIX_TIMESTAMP(CURRENT_TIMESTAMP()) where p_prepid = {}".format(status,self.pseudo_prepid)
                cur.execute(wtf)
            except MySQLdb.Error, e:
                raise Exception("ERROR %d IN CONNECTION: %s" % (e.args[0], e.args[1]))
            db.commit()
        finally:
            if db.open:
                db.close()

    def log_jid(self):
        jf="{}/.worker_{}.txt".format( self.scratch_dir, self.__class__.__name__) 
        msg = "run_locally=\t{}\njid=\t{}\nuser=\t{}@{}\npid=\t{}\npsid=\t{}\ndt=\t{}".format( 
          self.run_locally, (os.getenv("JOB_ID") if ("JOB_ID" in os.environ) else 'NULL'), 
          getpass.getuser(), socket.gethostname(), os.getpid(), self.pipeline_step_id, 
          datetime.datetime.isoformat(datetime.datetime.now()) 
        )
        with open(jf,"w") as f:
            f.write(msg)

    def set_dsm_status(self, status=0):

        S=90000+(10*self.pipeline_step_id)+status if status!=10 else 10

        db = get_connection("seqdb")
        try:
            cur = db.cursor()
            # horrid, but i am too sick of all this instability
            try:
                print("updating dsm to {} for {}".format(S,self.pseudo_prepid))
                cur.execute("update dragen_sample_metadata set is_merged = {} where pseudo_prepid = {}".format(S,self.pseudo_prepid))
            except MySQLdb.Error, e:
                raise Exception("ERROR %d IN CONNECTION: %s" % (e.args[0], e.args[1]))
            db.commit()
        finally:
            if db.open:
                db.close()
                
class JavaPipelineTask(GATKFPipelineTask):
    """Handle java's potentially large memory requirements by for example
    specifying that if one slot is requested, we get 24 GB memory/slot, or if 4
    are requested, we get 6/slot
    """
    javajob = True
    def __init__(self, *args, **kwargs):
        super(JavaPipelineTask, self).__init__(*args, **kwargs)
        self.mem = int(ceil(self.max_mem / self.n_cpu))
        if self.sample_type == "GENOME_AS_FAKE_EXOME":
            self.intervals_param = self.format_string(
                "-L {roche_bed}")
        else:
            self.intervals_param = ""

class FileExists(luigi.ExternalTask):
    fn = luigi.Parameter(
        description="the expected path to the DRAGEN aligned BAM")

    def output(self):
        return luigi.LocalTarget(self.fn)


class RealignerTargetCreator(JavaPipelineTask):
    priority = 1 # run before other competing steps with no requirement
    n_cpu = int(os.getenv("DEBUG_SLOTS")) if "DEBUG_SLOTS" in os.environ else 6
    max_mem = 24

    def pre_shell_commands(self):

        self.set_dsm_status(0) # just don't care anymore about this shite?!?
        self.update_sample_status('Realign Bam')
        self.log_jid()

        self.commands = [self.format_string(
            "{java} -Xmx{mem}g -jar {gatk} -R {ref_genome} -T RealignerTargetCreator "
            "-I {scratch_bam} -o {interval_list} -known {Mills1000g} "
            "-known {dbSNP} -nt {n_cpu} {intervals_param} {silly_arg}")]

    def requires(self):
        # return ValidateBAM(bam=self.scratch_bam, check_counts=self.sample_type != "CUSTOM_CAPTURE")
        return self.clone(EntryChecks)

    # def output(self):
        # luigi.LocalTarget(self.interval_list) 

class IndelRealigner(JavaPipelineTask):
    n_cpu = int(os.getenv("DEBUG_SLOTS")) if "DEBUG_SLOTS" in os.environ else 6
    max_mem = 24
    def pre_shell_commands(self):

        self.log_jid()
        self.set_dsm_status(0) 

        self.commands = [self.format_string(
            "{java} -Xmx{mem}g -jar {gatk} -R {ref_genome} -T IndelRealigner "
            "-I {scratch_bam} -o {realn_bam} -targetIntervals {interval_list} "
            "-maxReads 10000000 -maxInMemory 450000 -known {Mills1000g} "
            "-known {dbSNP} {intervals_param} {silly_arg}")]

    def requires(self):
        return self.clone(RealignerTargetCreator)

    # def output(self):
        # luigi.LocalTarget(self.realn_bam) 

class BaseRecalibrator(JavaPipelineTask):
    n_cpu = int(os.getenv("DEBUG_SLOTS")) if "DEBUG_SLOTS" in os.environ else 6
    max_mem = 24
    def pre_shell_commands(self):

        self.update_sample_status('Recal Bam')
        self.log_jid()
        self.set_dsm_status(0) 

        self.commands = [self.format_string(
            "{java} -Xmx{mem}g -jar {gatk} -R {ref_genome} "
            "-T BaseRecalibrator -I {realn_bam} "
            "-o {recal_table} -knownSites {Mills1000g} "
            "-knownSites {dbSNP} -nct {n_cpu} {intervals_param} {silly_arg}")]

    def requires(self):
        return self.clone(IndelRealigner)

class PrintReads(JavaPipelineTask):
    n_cpu = int(os.getenv("DEBUG_SLOTS")) if "DEBUG_SLOTS" in os.environ else 6
    max_mem = 24
    def pre_shell_commands(self):

        self.log_jid()
        self.set_dsm_status(0) 

        # --disable_indel_quals are necessary to remove BI and BD tags in the bam file
        self.commands = [self.format_string(
            "{java} -Xmx{mem}g -jar {gatk} -R {ref_genome} -T PrintReads "
            "-I {realn_bam} --disable_indel_quals "
            "-BQSR {recal_table} -o {recal_bam} -nct {n_cpu} {intervals_param} {silly_arg}")]

    def _run_post_success(self):
        os.remove(self.realn_bam)

    def requires(self):
        return self.clone(BaseRecalibrator)

    # def output(self):
        # luigi.LocalTarget(self.recal_bam) 
    # luigi.LocalTarget(self.realn_bam), luigi.LocalTarget(self.recal_table) ] 

class HaplotypeCaller(JavaPipelineTask):
    priority = 1 # run before other steps needing the recalibrated BAM
    n_cpu = int(os.getenv("DEBUG_SLOTS")) if "DEBUG_SLOTS" in os.environ else 8
    max_mem = 24
    def pre_shell_commands(self):

        self.arseholes()

        self.update_sample_status('Variant Calling')
        self.log_jid()
        self.set_dsm_status(0) 

        self.commands = [self.format_string(
            "{java} -Xms16g -Xmx34g -jar {gatk} -R {ref_genome} -T HaplotypeCaller "
            # dsth: making this sensible to avoid the bounds condition in 3.6...
            # dsth: updated gatk 3.6 to 2016-08-27-g667f78b
            # "-maxAltAlleles 3 "
            "-L {interval} -I {recal_bam} -o {gvcf} "
            "-stand_call_conf 20 -stand_emit_conf 20 --emitRefConfidence GVCF "
            "-GQB 5 -GQB 15 -GQB 20 -GQB 60 --variant_index_type LINEAR "
            "--variant_index_parameter 128000 --dbsnp {dbSNP} -nct {n_cpu}")]

    def post_shell_commands(self):

        self.set_dsm_status(1) 

        if os.system("gzip -t {}".format(self.gvcf))!=0:
            self.set_dsm_status(3) 
            raise VCFCheckException('corrupted block compression file')

        chr22check = subprocess.Popen('/nfs/goldstein/software/tabix-0.2.6/tabix {} 22 | head -10 | wc -l'.format(self.gvcf),shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if chr22check.stdout.readline().rstrip() != '10': # f'ing lazy
            self.set_dsm_status(3) 
            raise VCFCheckException('truncated file')

        ##### this clearly isn't very exhaustive...
        errors = check_vcf(self.gvcf)
        if errors:
            self.set_dsm_status(3) 
            raise VCFCheckException("\n".join(errors))

        self.set_dsm_status(2) 

    def requires(self):
        return self.clone(PrintReads)

class GenotypeGVCFs(GATKFPipelineTask):
    priority = 1
    def pre_shell_commands(self):

        self.arseholes()

        self.update_sample_status('Genotyping')
        self.log_jid()
        self.set_dsm_status(0) 

        self.commands = [self.format_string(
            "{java} -jar {gatk} -R {ref_genome} -T GenotypeGVCFs "
            "-L {interval} -o {vcf} -stand_call_conf 20 -stand_emit_conf 20 -V {gvcf}")]

    def post_shell_commands(self):

        self.set_dsm_status(1) 

        errors = check_vcf(self.vcf, self.check_variant_counts)
        if errors:
            self.set_dsm_status(3) 
            msg="{} QC Errors ({}) : Last Error {}".format(len(errors),self.pipeline_step_id,errors[-1])
            self.set_dsm_status(10000) 
            self.update_sample_status(msg)
            exit(0)
            raise QCError(msg)
        self.set_dsm_status(2) 

    def requires(self):
        return self.clone(HaplotypeCaller)

class SelectVariantsSNP(GATKFPipelineTask):
    priority = 1
    def pre_shell_commands(self):

        self.arseholes()

        self.log_jid()
        self.set_dsm_status(0) 

        self.commands = [self.format_string(
            "{java} -jar {gatk} -R {ref_genome} -T SelectVariants "
            "-L {interval} -V {vcf}  -selectType SNP -o {snp_vcf}")]

    def post_shell_commands(self):

        self.set_dsm_status(1) 

        errors = check_vcf(self.snp_vcf)
        if errors:
            self.set_dsm_status(3) 
            raise VCFCheckException("\n".join(errors))
        self.set_dsm_status(2) 

    def requires(self):
        return self.clone(GenotypeGVCFs)

class SelectVariantsINDEL(GATKFPipelineTask):
    def pre_shell_commands(self):

        self.log_jid()
        self.set_dsm_status(0) 

        self.commands = [self.format_string(
            "{java} -jar {gatk} -R {ref_genome} -T SelectVariants "
            "-L {interval} -V {vcf}  -selectType INDEL -o {indel_vcf}")]

    def post_shell_commands(self):

        self.set_dsm_status(1) 

        errors = check_vcf(self.indel_vcf)
        if errors:
            self.set_dsm_status(3) 
            raise VCFCheckException("\n".join(errors))
        self.set_dsm_status(2) 

    def requires(self):
        return self.clone(GenotypeGVCFs)

class VariantRecalibratorSNP(JavaPipelineTask):
    max_mem = 16
    def pre_shell_commands(self):
        ## Running both Exomes and Genomes the same way, i.e. we will exlcude DP, omni is set to be a truth set, also added in ExAC SNPs as a training
        ## set which can contain both TP and FP with the same prior likelihood as the 1000G training set. 
        ## See these links : 
        ## For exomes : 
        ## https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php
        ## For genomes : 
        ## https://software.broadinstitute.org/gatk/guide/article?id=1259

        self.log_jid()
        self.set_dsm_status(0) 

        self.commands = [self.format_string(
            "{java} -Xmx{mem}g -jar {gatk} -R {ref_genome} -T VariantRecalibrator "
            "-L {interval}  --input {snp_vcf} -an QD -an FS -an SOR -an MQ "
            "-an MQRankSum -an ReadPosRankSum "
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

        self.arseholes()

        self.log_jid()
        self.set_dsm_status(0) 

        self.commands = [self.format_string(
            '{java} -jar {gatk} -R {ref_genome} -T VariantFiltration '
            '-L {interval} -V {snp_vcf} --filterExpression '
            '"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || '
            'ReadPosRankSum < -8.0" --filterName "SNP_filter" -o {snp_filtered}')]

    def post_shell_commands(self):

        self.set_dsm_status(1) 

        errors = check_vcf(self.snp_filtered)
        if errors:
            self.set_dsm_status(3) 
            raise VCFCheckException("\n".join(errors))
        self.set_dsm_status(2) 

    def requires(self):
        return self.clone(SelectVariantsSNP)

class VariantRecalibratorINDEL(JavaPipelineTask):
    max_mem = 16
    def pre_shell_commands(self):

        self.arseholes()

        self.log_jid()
        self.set_dsm_status(0) 

        self.commands = [self.format_string(
            "{java} -Xmx{mem}g -jar {gatk} -R {ref_genome} "
            "-T VariantRecalibrator -L {interval} "
            "--input {indel_vcf} -an QD -an FS -an SOR -an MQRankSum "
            "-an ReadPosRankSum -mode INDEL --maxGaussians 4 "
            "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 "
            "-recalFile {indel_recal} -tranchesFile {indel_tranches} "
            "-resource:mills,known=true,training=true,truth=true,prior=12.0 {Mills1000g}")]

    def requires(self):
        return self.clone(SelectVariantsINDEL)

class ApplyRecalibrationSNP(GATKFPipelineTask):
    def pre_shell_commands(self):

        self.arseholes()

        self.log_jid()
        self.set_dsm_status(0) 

        self.commands = [self.format_string(
            "{java} -jar {gatk} -R {ref_genome} -T ApplyRecalibration "
            "-L {interval}  -input {snp_vcf} "
            "-tranchesFile {snp_tranches} -recalFile {snp_recal} "
            "-o {snp_filtered} --ts_filter_level 90.0 -mode SNP")]

    def post_shell_commands(self):

        self.set_dsm_status(1) 

        errors = check_vcf(self.snp_filtered)
        if errors:
            self.set_dsm_status(3) 
            raise VCFCheckException("\n".join(errors))
        self.set_dsm_status(2) 

    def requires(self):
        return self.clone(VariantRecalibratorSNP)

class ApplyRecalibrationINDEL(GATKFPipelineTask):
    def pre_shell_commands(self):

        self.arseholes()

        self.log_jid()
        self.set_dsm_status(0) 

        self.commands = [self.format_string(
            "{java} -jar {gatk} -R {ref_genome} -T ApplyRecalibration "
            "-L {interval}  -input {indel_vcf} "
            "-tranchesFile {indel_tranches} -recalFile {indel_recal} "
            "-o {indel_filtered} --ts_filter_level 90.0 -mode INDEL")]

    def post_shell_commands(self):

        self.set_dsm_status(1) 

        errors = check_vcf(self.indel_filtered)
        if errors:
            self.set_dsm_status(3) 
            raise VCFCheckException("\n".join(errors))
        self.set_dsm_status(2) 

    def requires(self):
      return self.clone(VariantRecalibratorINDEL)

class VariantFiltrationINDEL(GATKFPipelineTask):
    def pre_shell_commands(self):

        self.arseholes()

        self.log_jid()
        self.set_dsm_status(0) 

        self.commands = [self.format_string(
            '{java} -jar {gatk} -R {ref_genome} -T VariantFiltration '
            '-L {interval} -V {indel_vcf} --filterExpression '
            '"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName '
            '"INDEL_filter"  -o {indel_filtered}')]

    def post_shell_commands(self):

        self.set_dsm_status(1) 

        errors = check_vcf(self.indel_filtered)
        if errors:
            self.set_dsm_status(3) 
            raise VCFCheckException("\n".join(errors))
        self.set_dsm_status(2) 

    def requires(self):
      return self.clone(SelectVariantsINDEL)

class CombineVariants(GATKFPipelineTask):
    def pre_shell_commands(self):

        self.arseholes()

        self.log_jid()
        self.set_dsm_status(0) 

        """Merges SNP and INDEL vcfs.  Using the filtered SNP vcf header as a
        base, variant type/sample type specific ##FILTERs are added to the header.
        After finishing reading through the header the SNP vcf is read again with 
        only variants being outputted.  The same happens with the INDEL vcf.
        The resulting vcf is then sorted producing the analysisReady.vcf.
        """
        filter_encountered = False
        info_encountered = False
        
        with open(self.tmp_vcf, "w") as vcf_out:
            for line in open(self.snp_filtered).readlines():
                if line.startswith("#"):
                    if line.startswith("##FILTER") and not filter_encountered:
                        filter_encountered = True
                        with open(self.indel_filtered) as indel_header:
                            prefix = ("##FILTER=<ID=VQSRTrancheINDEL" if self.sample_type == "GENOME" else "##FILTER=<ID=INDEL_filter,")
                            indel_lines = []

                            for indel_line in indel_header:
                                if indel_line.startswith(prefix):
                                    indel_lines.append(indel_line)
                                elif not indel_line.startswith("#"):
                                    break

                            if not indel_lines:
                                self.set_dsm_status(3) 
                                raise ValueError(
                                    "Could not parse indel filter/VQSR line(s) "
                                    "from filtered indel VCF")
                            for indel_line in indel_lines:
                                vcf_out.write(indel_line)

                    if line.startswith("##INFO") and not info_encountered:
                        info_encountered = True
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
            "{java} -XX:ParallelGCThreads=1 -jar {picard} SortVcf I={tmp_vcf} O={final_vcf}"))
        self.commands.append(self.format_string("{bgzip} -f {final_vcf}"))
        self.commands.append(self.format_string(
            "{tabix} -f {final_vcf_gz}"))

    def post_shell_commands(self):

        self.set_dsm_status(1) 

        errors = check_vcf(self.final_vcf_gz, self.check_variant_counts)
        if errors:
            self.set_dsm_status(3) 
            raise VCFCheckException("\n".join(errors))
        self.set_dsm_status(2) 

    def requires(self):
        if self.sample_type in ("EXOME", "GENOME", "GENOME_AS_FAKE_EXOME"):
            yield self.clone(ApplyRecalibrationSNP)
            if self.sample_type in ("EXOME", "GENOME_AS_FAKE_EXOME"):
                yield self.clone(VariantFiltrationINDEL)
            else:
                yield self.clone(ApplyRecalibrationINDEL)
        elif self.sample_type == "CUSTOM_CAPTURE":
            yield self.clone(VariantFiltrationSNP)
            yield self.clone(VariantFiltrationINDEL)
        else:
            self.set_dsm_status(3) 
            raise ValueError(
                "Sample type: {} not supported".format(self.sample_type))

class RBP(GATKFPipelineTask):
    priority = 1
    def pre_shell_commands(self):

        self.arseholes()

        self.log_jid()
        self.set_dsm_status(0) 

        self.commands = [self.format_string(
            "{java} -XX:ParallelGCThreads=1 -jar {gatk} -T ReadBackedPhasing -R {ref_genome} "
            "-I {recal_bam} --variant {final_vcf_gz} -o {phased_vcf} "
            "--phaseQualityThresh 20.0 --maxGenomicDistanceForMNP 2 "
            "--enableMergePhasedSegregatingPolymorphismsToMNP "
            "-U ALLOW_SEQ_DICT_INCOMPATIBILITY")]

    def post_shell_commands(self):

        self.set_dsm_status(1) 

        errors = check_vcf(self.phased_vcf, self.check_variant_counts)
        if errors:
            self.set_dsm_status(3) 
            raise VCFCheckException("\n".join(errors))
        self.set_dsm_status(2) 

    def requires(self):
        return self.clone(CombineVariants), self.clone(PrintReads)

class FixMergedMNPInfo(GATKFPipelineTask):
    """ The MNPs phased by RBP are missing the DP,AD and other annotation info
    this task which check variants which are missing these information, fetch
    the value for corresponding to that variant site from the vcf prior to RBP
    and add these information to the fixed vcf """
    def post_shell_commands(self):

        self.set_dsm_status(1) 

        self.fix_phased_vcf()
        errors = check_vcf(self.fixed_vcf, self.check_variant_counts)
        if errors:
            self.set_dsm_status(3) 
            raise VCFCheckException("\n".join(errors))
        self.set_dsm_status(2) 

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

        self.arseholes()

        self.update_sample_status('Annotating Variants')
        self.log_jid()
        self.set_dsm_status(0) 

        self.shell_options.update(
            {"stdout":None, "stderr":self.log_file, "shell":True})
        self.commands = [self.format_string(
            "/nfs/goldstein/software/jdk1.8.0_05/bin/java "
            "-XX:ParallelGCThreads=1 -Xmx6g -jar {clineff} "
            "-c {clineff_cfg} -v -db {annotatedbSNP} GRCh37.87 {fixed_vcf} "
            "> {annotated_vcf}"),
            # bc, 10/11/17: ClinEff has been observed to silently fail, i.e.
            # produce an emtpy VCF, so check for that
            self.format_string("if [ -s {annotated_vcf} ]; then {bgzip} -f {annotated_vcf}; "
                               "else echo 1>&2 {annotated_vcf} is empty!; exit 1; fi"),
            #self.format_string("{bgzip} -f {annotated_vcf}"),
            self.format_string("{tabix} -f {annotated_vcf_gz}")]

    def post_shell_commands(self):

        self.set_dsm_status(1) 

        errors = check_vcf(self.annotated_vcf_gz, self.check_variant_counts)
        if errors:
            self.set_dsm_status(3) 
            raise VCFCheckException("\n".join(errors))
        self.set_dsm_status(2) 

    def requires(self):
        return self.clone(FixMergedMNPInfo)

class SubsetVCF(GATKFPipelineTask):
    """ For SRR/Any future deemed to be low quality sequenced samples, subset the vcf to only the capturekit regions """
    def pre_shell_commands(self):

        self.log_jid()
        self.set_dsm_status(0) 

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

    def post_shell_commands(self):

        self.set_dsm_status(1) 

        errors = check_vcf(self.annotated_vcf_gz, self.check_variant_counts)
        if errors:
            self.set_dsm_status(3) 
            raise VCFCheckException("\n".join(errors))
        self.set_dsm_status(2) 

    def requires(self):
        return self.clone(AnnotateVCF)

####### there's some pretty strange ordering here?!?
class PreArchiveChecks(GATKFPipelineTask):

    n_cpu = int(os.getenv("DEBUG_SLOTS")) if "DEBUG_SLOTS" in os.environ else 2
    prept_start_message = "PreArchiveChecks"
    prept_completed_message = "PreArchiveChecks"
    # prept_start_message = prept_completed_message = __class__.__name__

    dont_remove_tmp_dir = False # remove the temporary directory iff this task succeeds
    dont_remove_tmp_dir_if_failure = True # don't remove if it fails

    def requires(self):
        if self.sample_name.upper().startswith('SRR'):
            yield self.clone(SubsetVCF)
        # yield self.clone(SplitAndSubsetDPBins)
        yield self.clone(UpdateSeqdbMetrics)

    def __init__(self, *args, **kwargs):

        super(PreArchiveChecks, self).__init__(*args, **kwargs)

        if (self.sample_name.upper().startswith('PGMCLIN') or self.sample_name.upper().startswith('PGMVIP')):
            self.base_dir = os.path.join( "/nfs/pgmclin/ALIGNMENT/BUILD37/DRAGEN/", self.sample_type, self.name_prep )
        else:
            self.base_dir = os.path.join( self.config_parameters["base"], self.sample_type, self.name_prep )

        self.check_counts = self.sample_type != "CUSTOM_CAPTURE"

    def pre_shell_commands(self):

        self.log_jid()
        self.set_dsm_status(0) 

        self.script_dir = os.path.join(self.scratch_dir, "scripts")

        self.pipeline_tarball   = ( "{scratch_dir}/{name_prep}.pipeline_data.tar.gz".format( scratch_dir=self.scratch_dir, name_prep=self.name_prep) )
        self.cvg_tarball        = ( "{cov_dir}/coverage.tar.gz".format(cov_dir=self.cov_dir) )
        self.gq_tarball         = ( "{gq_dir}/gq.tar.gz".format(gq_dir=self.gq_dir) )
        self.raw_coverage       = os.path.join( self.scratch_dir, "{}.coverage_bins".format(self.name_prep) )

        ####### DISABLE_ME
        # if ( datetime.datetime.now().strftime('%Y-%m-%d')=='2018-06-01' or datetime.datetime.now().strftime('%Y-%m-%d')=='2014-01-02' ) and os.path.isdir(self.base_dir): 
        print(datetime.datetime.now().strftime('%Y-%m-%d'))

        if ( os.path.isdir(self.base_dir) ): 
        # if ( datetime.datetime.now().strftime('%Y-%m-%d')=='2018-07-08') and os.path.isdir(self.base_dir): 
            print("disable this : so bored of the messes - wiping for now but really ought to rsync again at min but really check WTF is going on!?!?")
            cmd = "echo '{}' >> /home/dh2880/arseholes_wipe_these.txt".format(self.base_dir)
            os.system(cmd)
            try:
                print("attempting to remove '{}'".format(self.base_dir))
                rmtree(self.base_dir)
            except:
                print("unable to clean up archive dir : {} ".format(self.base_dir))
                self.set_dsm_status(10001) 
                exit(0)

        if os.path.isdir(self.base_dir): 
            ### hacky but just don't care anymore!?!
            print("archive dir exists : {} ".format(self.base_dir))
            self.set_dsm_status(10000) 
            exit(0)
            raise ArchiveDirectoryAlreadyExists( "the archive location, '{}', already exists".format(self.base_dir) )

        frameinfo = getframeinfo(currentframe())
        vcf_errors = check_vcf(self.annotated_vcf_gz, self.check_variant_counts)
        if vcf_errors:
            if "DEBUG_INTERVALS" not in os.environ:
                self.set_dsm_status(3) 
                raise VCFCheckException( "\n".join(vcf_errors) )

        frameinfo = getframeinfo(currentframe())
        # try:
        bam_errors = check_bam( self.recal_bam, self.check_counts )
        if bam_errors:
            if "DEBUG_INTERVALS" not in os.environ:
                self.set_dsm_status(3) 
                raise BAMCheckException( "\n".join(bam_errors) )

        frameinfo = getframeinfo(currentframe())

        if os.path.exists(self.pipeline_tarball):
            # raise ValueError("\n\ntar file already exists '{}'".format(self.pipeline_tarball))
            try:
                os.remove(self.pipeline_tarball)
            except:
                raise ValueError("\n\nunable to remove old tar file '{}'".format(self.pipeline_tarball))

        with tarfile.open(  self.pipeline_tarball, "w:gz"   ) as tar:
            for d in (self.script_dir, self.log_dir):
                tar.add( d, arcname=os.path.basename(d) )
            for txt in glob("{scratch_dir}/*.txt".format(scratch_dir=self.scratch_dir)):
                tar.add( txt, arcname=os.path.basename(txt) )

        with tarfile.open(  self.cvg_tarball, "w:gz"        ) as tar:
            for cvg_file in glob( "{cov_dir}/*.txt".format(cov_dir=self.cov_dir)):
                tar.add( cvg_file, arcname=os.path.basename(cvg_file) )

        with tarfile.open(  self.gq_tarball, "w:gz"         ) as tar:
            for gq_file in glob( "{gq_dir}/*.txt".format(gq_dir=self.gq_dir)):
                tar.add( gq_file, arcname=os.path.basename(gq_file) )

        frameinfo = getframeinfo(currentframe())

    def run(self):
        try:
            super(PreArchiveChecks, self).run()
        except Exception, e:
            if (type(e) is not ArchiveDirectoryAlreadyExists
                and os.path.isdir(self.base_dir)):
                # clean up the directory if it was created and an error was
                # generated
                try:
                    rmtree(self.base_dir)
                except:
                    print("unable to clean up archive dir : {} ".format(self.base_dir))
                    self.set_dsm_status(10001) 
                    exit(0)

class ArchiveSample(GATKFPipelineTask):

    n_cpu = int(os.getenv("DEBUG_SLOTS")) if "DEBUG_SLOTS" in os.environ else 2
    prept_start_message = "ArchiveSample"
    prept_completed_message = "ArchiveSample" 

    dont_remove_tmp_dir = False # remove the temporary directory iff this task succeeds
    dont_remove_tmp_dir_if_failure = True # don't remove if it fails

    def requires(self):
        yield self.clone(PreArchiveChecks)

    def __init__(self, *args, **kwargs):
        super(ArchiveSample, self).__init__(*args, **kwargs)
        if (self.sample_name.upper().startswith('PGMCLIN') or
            self.sample_name.upper().startswith('PGMVIP')):
            self.base_dir = os.path.join(
                "/nfs/pgmclin/ALIGNMENT/BUILD37/DRAGEN/",
                self.sample_type, self.name_prep)
        else:
            self.base_dir = os.path.join(
                self.config_parameters["base"], self.sample_type, self.name_prep)
        self.check_counts = self.sample_type != "CUSTOM_CAPTURE"

    def pre_shell_commands(self):

        self.update_sample_status('Archiving')
        self.log_jid()
        self.set_dsm_status(0) 

        self.script_dir = os.path.join(self.scratch_dir, "scripts")

        self.pipeline_tarball   = ( "{scratch_dir}/{name_prep}.pipeline_data.tar.gz".format( scratch_dir=self.scratch_dir, name_prep=self.name_prep) )
        self.cvg_tarball        = ( "{cov_dir}/coverage.tar.gz".format(cov_dir=self.cov_dir) )
        self.gq_tarball         = ( "{gq_dir}/gq.tar.gz".format(gq_dir=self.gq_dir) )
        self.raw_coverage       = os.path.join( self.scratch_dir, "{}.coverage_bins".format(self.name_prep) )

        if not os.path.isdir(self.base_dir): 

            frameinfo = getframeinfo(currentframe())

            self.archive_helper() # just get the proper list

            if not os.path.isdir(self.base_dir): 
                os.makedirs(self.base_dir)

            for data_file in self.data_to_copy:
                f = self.format_string("{}".format(data_file))
                ### WTF : suffix changed for some dragen initial bams - merging?!?
                if data_file=="{scratch_bam}.bai" and not os.path.isfile(f):
                    n=f.replace('bam.bai','bai')
                    os.rename(n,f)
                if not os.path.isfile(f):
                    raise ValueError("\n\nfile {} does not exist\n".format(f))
                # print(self.format_string("copy {} to {}".format(data_file,self.base_dir)))
                # print("copy {} to {}".format(data_file,self.base_dir))
                self.commands.append(self.format_string( "rsync -grlt --inplace --partial " + data_file + " {base_dir}") )

            # re-check the archived data to ensure its integrity before deleting the
            # scratch directory
            frameinfo = getframeinfo(currentframe())
            print('will copy', frameinfo.filename, frameinfo.lineno)
        else:
            frameinfo = getframeinfo(currentframe())
            self.set_dsm_status(10000) 
            print("this shouldn't exist '{}'".format(self.base_dir))
            exit(0)
            # raise ValueError('\n\nthe archive dir is present {}'.format(self.base_dir))

    def post_shell_commands(self):

        self.set_dsm_status(1) 

        for data_file in self.data_to_copy:
            scratch_fn = self.format_string(data_file)
            archived_fn = os.path.join(self.base_dir, os.path.basename(scratch_fn))
            if os.path.getsize(scratch_fn) != os.path.getsize(archived_fn):
                self.set_dsm_status(3) 
                raise RsyncException("Data size of original file ({}) does not match the archived version ({})!".format( scratch_fn, archived_fn) )

##################################################################
# this really shouldn't be here but i just don't want the bother...
##################################################################

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

        vcf_errors = check_vcf(
            os.path.join(self.base_dir, os.path.basename(self.annotated_vcf_gz)),
            self.check_variant_counts)
        if vcf_errors:
            if "DEBUG_INTERVALS" not in os.environ:
                raise VCFCheckException("\n".join(vcf_errors))
        bam_errors = check_bam(
            os.path.join(self.base_dir, os.path.basename(self.recal_bam)),
            self.check_counts)
        if bam_errors:
            if "DEBUG_INTERVALS" not in os.environ:
                raise BAMCheckException("\n".join(bam_errors))
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

        self.set_dsm_status(10)

    def run(self):
        try:
            super(ArchiveSample, self).run()
        except Exception, e:
            if (type(e) is not ArchiveDirectoryAlreadyExists
                and os.path.isdir(self.base_dir)):
                # clean up the directory if it was created and an error was
                # generated
                print("something went wrong so cleaning up archive dir '{}'".format(self.base_dir))
                try:
                    rmtree(self.base_dir)
                except:
                    print("unable to clean up archive dir : {} ".format(self.base_dir))
                    self.set_dsm_status(10001) 
                    exit(0)
            raise

class CoverageBinning(GATKFPipelineTask):
    """Call Dan's program to bin coverage where we require some minimum
    standards in base quality and mapping quality, and only output records for
    blocks which have at least one base meeting the specified binned depth value
    """
    bin_program = luigi.InputFileParameter(
        description="The path to the binning program to run")
    mmq = luigi.IntParameter(
        description="Require a minimum base quality score of this value for a "
        "base in a read to be counted")
    mmb = luigi.IntParameter(
        description="Require a minimum mapping quality score of this value for "
        "a given read to be counted")
    mbd = luigi.ChoiceParameter(
        choices=["a", "b", "c", "d", "e", "f", "g"],
        description="Require a block to have at least one base of the given "
        "binned coverage value to be output")

    def pre_shell_commands(self):

        self.set_dsm_status(0) 

        self.shell_options.update(
            {"stdout":os.path.join(
                self.scratch_dir, "{}.coverage_bins".format(self.name_prep)),
             "stderr":self.log_file})
        self.commands = [self.format_string(
            "{bin_program} {recal_bam} {mmq} {mmb} {mbd}")]

        print("scratch_dir='{}'".format(self.name_prep))

    def post_shell_commands(self):
        self.set_dsm_status(1) 
        of=os.path.join(self.scratch_dir,"{}.coverage_bins".format(self.name_prep))
        if os.system('perl -ne "exit 1 if(\$_=~/[\\x00-\\x08\\x7F-\\xFF]/)" {}'.format(of))!=0:
            self.set_dsm_status(3) 
            raise ValueError("\n\nthere appears to be an issue with the output file {}\n".format(of))

        p = subprocess.Popen('grep -c -P "[\\x00-\\x08\\x7F-\\xFF]" {}'.format(of), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        probs=p.stdout.readline().rstrip()
        if probs != "0":
            raise ValueError("\n\nthere appears to be an issue with the output file {}\n".format(of))

        self.set_dsm_status(2) 

    def requires(self):
        return self.clone(PrintReads)

class SplitAndSubsetDPBins(GATKFPipelineTask):
    """Split the DP bins by chromosome and subset them to the blocks of interest
    if they aren't genomes
    """
    def pre_shell_commands(self):

        self.set_dsm_status(0) 

        if not os.path.isdir(self.cov_dir):
            os.makedirs(self.cov_dir)
        if self.sample_type == "GENOME":
            dp_blocks_fn = None
        else:
            seqdb = get_connection("seqdb")
            try:
                seq_cur = seqdb.cursor()
                seq_cur.execute(GET_DP_BLOCKS_FILE.format(
                    capture_kit=self.capture_kit))
                row = seq_cur.fetchone()
                if row:
                    dp_blocks_fn = row[0]
                else:
                    self.set_dsm_status(3) 
                    raise ValueError(
                        "Could not find DP blocks file name for capture kit {}".
                        format(self.capture_kit))
                if not os.path.isfile(dp_blocks_fn):
                    self.set_dsm_status(3) 
                    raise OSError("DP blocks file does not exist for capture kit {}".
                                  format(self.capture_kit))
            finally:
                if seqdb.open:
                    seqdb.close()

        coverage_files = dict([chrom, os.path.join(self.cov_dir,
                          "{name_prep}_coverage_binned_1000_chr{chrom}.txt".format(
                              name_prep=self.name_prep, chrom=chrom))]
                              for chrom in CHROMs)
        for coverage_fn in coverage_files.itervalues():
            with open(coverage_fn, "w"): pass
        if dp_blocks_fn:
            blocks_to_retain = defaultdict(set)
            with open(dp_blocks_fn) as dp_blocks_fh:
                for line in dp_blocks_fh:
                    chrom, block_id = line.strip().split(":")
                    blocks_to_retain[chrom].add(block_id)
        coverage_out_fh = None
        prev_chrom = None
        with open(os.path.join( self.scratch_dir, "{}.coverage_bins".format(self.name_prep)) ) as coverage_fh:
            for x, line in enumerate(coverage_fh):
                chrom, block_id = line.split("\t")[:2]
                if dp_blocks_fn:
                    output_record = block_id in blocks_to_retain[chrom]
                    if block_id in blocks_to_retain[chrom] or chrom == "MT":
                        output_record = True
                        if prev_chrom != chrom:
                            coverage_out_fh = open(coverage_files[chrom], "w")
                            prev_chrom = chrom
                    else:
                        output_record = False
                else:
                    output_record = chrom in CHROMs
                if output_record:
                    if prev_chrom != chrom:
                        if coverage_out_fh:
                            coverage_out_fh.close()
                        coverage_out_fh = open(coverage_files[chrom], "w")
                        prev_chrom = chrom
                    coverage_out_fh.write(line)
        if coverage_out_fh:
            coverage_out_fh.close()

    def requires(self):
        return self.clone(CoverageBinning)
   
class GQBinning(GATKFPipelineTask):
    def pre_shell_commands(self):

        self.set_dsm_status(0) 

        try:
            if not os.path.isdir(self.gq_dir):
                os.makedirs(self.gq_dir)
        except OSError,e:
            if e.errno != 17:
                self.set_dsm_status(3) 
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

        self.set_dsm_status(0) 

        self.alignment_metrics_raw = os.path.join(
            self.scratch_dir,
            "{sample_name}.{pseudo_prepid}.alignment.metrics.raw".format(
                sample_name=self.sample_name, pseudo_prepid=self.pseudo_prepid))
        cmd = self.format_string(
            "{java} -XX:ParallelGCThreads=1 -jar {picard} CollectAlignmentSummaryMetrics "
            "TMP_DIR={scratch_dir} "
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

        self.set_dsm_status(0) 

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
        ## An intermediate parsed file is created for extracting info for db update later, this remains the same for exomes and genomes
        if self.sample_type == "GENOME":
            ## Define shell commands to be run
            cvg_cmd = ("{java} -XX:ParallelGCThreads=1 -jar "
                       "{picard} CollectWgsMetrics VALIDATION_STRINGENCY=LENIENT "
                       "R={ref_genome} I={recal_bam} INTERVALS={file_target} "
                       "O={raw_output_file} MQ=20 Q=10 >> {log_file} 2>&1")
            ## Run on the ccds regions only
            cvg_cmd1 = self.format_string(
                cvg_cmd, raw_output_file=self.raw_output_file_ccds,
                file_target=self.config_parameters["target_file"])
            ## Run across the genome
            cvg_cmd2 = self.format_string(
                "{java} -XX:ParallelGCThreads=1 -jar "
                "{picard} CollectWgsMetrics R={ref_genome} I={recal_bam} "
                "O={raw_output_file} MQ=20 Q=10 >> {log_file} 2>&1")
            ## Run on X and Y Chromosomes only (across all regions there not just ccds)
            cvg_cmd3 = self.format_string(cvg_cmd,
                file_target=self.config_parameters["target_file_X"],
                raw_output_file=self.raw_output_file_X)
            cvg_cmd4 = self.format_string(cvg_cmd,
                file_target=self.config_parameters["target_file_Y"],
                raw_output_file=self.raw_output_file_Y)
        else:
            db = get_connection("seqdb")
            try:
                cur = db.cursor()
                query = ("SELECT region_file_lsrc FROM captureKit "
                         "INNER JOIN dragen_sample_metadata ON prepT_name=dragen_sample_metadata.capture_kit "
                         "WHERE pseudo_prepid = {0} AND chr = 'X'".format(self.pseudo_prepid))
                cur.execute(query)
                row = cur.fetchone()
                if row:
                    self.capture_file_X = row[0]
                elif self.sample_type == "CUSTOM_CAPTURE":
                    # we don't require X coverage for these
                    self.capture_file_X = None
                else:
                    self.set_dsm_status(3) 
                    raise ValueError("Couldn't find BED file for coverage on "
                                    "X chromosome")
                query = ("SELECT region_file_lsrc FROM captureKit "
                         "INNER JOIN dragen_sample_metadata ON prepT_name=dragen_sample_metadata.capture_kit "
                         "WHERE pseudo_prepid = {0} AND chr = 'Y'".format(self.pseudo_prepid))
                cur.execute(query)
                row = cur.fetchone()
                if row:
                    self.capture_file_Y = row[0]
                elif self.sample_type == "CUSTOM_CAPTURE":
                    self.capture_file_Y = None
                else:
                    self.set_dsm_status(3) 
                    raise ValueError("Couldn't find BED file for coverage on "
                                     "Y chromosome")
            finally:
                if db.open:
                    db.close()

            self.calc_capture_specificity = 1 ## We need to run an additional picard module for calculating capture specificity
            cvg_cmd = ("{java} -XX:ParallelGCThreads=1 -jar "
                       "{gatk} -T DepthOfCoverage -mbq 10 -mmq 20 "
                       "--omitIntervalStatistics --omitDepthOutputAtEachBase "
                       "--omitLocusTable -R {ref_genome} -I {recal_bam} "
                       "-ct 5 -ct 10 -ct 15 -ct 20 -L {file_target} "
                       "-o {raw_output_file} >> {log_file} 2>&1")
            ## Run across ccds regions
            cvg_cmd1 = self.format_string(cvg_cmd,
                file_target=self.config_parameters["target_file"],
                raw_output_file=self.raw_output_file_ccds)
            ## Run across capture kit regions
            cvg_cmd2 = self.format_string(cvg_cmd,
                file_target=self.capture_kit_bed,
                raw_output_file=self.raw_output_file)
            ## Run on X and Y Chromosomes only (across all regions there not just ccds)
            if self.capture_file_X:
                cvg_cmd3 = self.format_string(
                    cvg_cmd, file_target=self.capture_file_X,
                    raw_output_file=self.raw_output_file_X)
            else:
                cvg_cmd3 = None
            if self.capture_file_Y:
                cvg_cmd4 = self.format_string(
                    cvg_cmd, file_target=self.capture_file_Y,
                    raw_output_file=self.raw_output_file_Y)
            else:
                cvg_cmd4 = None
            ## Run PicardHsMetrics for Capture Specificity
            self.output_bait_file = os.path.join(self.scratch_dir, "targets.interval")
            cvg_cmd5 = self.format_string(
                "{java} -XX:ParallelGCThreads=1 -jar {picard} "
                "CollectHsMetrics BI={output_bait_file} TI={output_bait_file} "
                "VALIDATION_STRINGENCY=SILENT "
                "METRIC_ACCUMULATION_LEVEL=ALL_READS "
                "I={recal_bam} O={raw_output_file_cs} MQ=20 Q=10 >> {log_file} 2>&1")
            ## DepthOfCoverage output has an extra suffix at the end
            self.raw_output_file_ccds = self.raw_output_file_ccds + '.sample_summary'
            self.raw_output_file = self.raw_output_file + '.sample_summary'
            self.raw_output_file_X = self.raw_output_file_X + '.sample_summary'
            self.raw_output_file_Y = self.raw_output_file_Y + '.sample_summary'

        parser_cmd = "grep -v '^#' {raw_output_file} | awk -f {transpose_awk} > {output_file}"
        parser_cmd1 = self.format_string(
            parser_cmd, raw_output_file=self.raw_output_file_ccds, output_file=self.output_file_ccds)
        parser_cmd2 = self.format_string(
            parser_cmd, raw_output_file=self.raw_output_file, output_file=self.output_file)
        self.commands = [cvg_cmd1, parser_cmd1, cvg_cmd2, parser_cmd2]
        if cvg_cmd3:
            self.commands.append(cvg_cmd3)
            self.commands.append(self.format_string(
                parser_cmd, raw_output_file=self.raw_output_file_X,
                output_file=self.output_file_X))
        if cvg_cmd4:
            self.commands.append(cvg_cmd4)
            self.commands.append(self.format_string(
                parser_cmd, raw_output_file=self.raw_output_file_Y,
                output_file=self.output_file_Y))
        if self.sample_type != "GENOME": ## i.e. Exome or CustomCapture
            self.commands.append(self.format_string(
                "{java} -XX:ParallelGCThreads=1 -jar {picard} BedToIntervalList I={capture_kit_bed} "
                "SD={seqdict_file} OUTPUT={output_bait_file} >> {log_file} 2>&1"))
            self.commands.append(cvg_cmd5)
            self.commands.append(self.format_string(
                parser_cmd, raw_output_file=self.raw_output_file_cs,
                output_file=self.output_file_cs))
        self.shell_options.update({"stdout":None, "stderr":None, "shell":True})

    def requires(self):
        yield self.clone(PrintReads)

class EntryChecks(GATKFPipelineTask):

    """ Parse Duplicate Metrics from dragen logs """

    def __init__(self, *args, **kwargs):

        super(EntryChecks, self).__init__(*args, **kwargs)
        self.run_locally=True
        self.dragen_log = os.path.join( self.log_dir, self.name_prep + ".dragen.out" )
        self.picard_log = os.path.join( self.scratch_dir, self.name_prep + ".metrics_duplication.txt" )
        self.duplicates_file = os.path.join( self.scratch_dir, self.name_prep + ".duplicates.txt" )
        self.merge_intermdiate = os.path.join( self.scratch_dir, self.name_prep + ".merge.bam" )

    def pre_shell_commands(self):

        self.log_jid()
        self.set_dsm_status(0) 

        if not os.path.exists(self.scratch_bam):
            self.set_dsm_status(3) 
            raise ValueError("There's no bam file ({})".format(self.scratch_bam))
            # with self.output().open("w"): pass

        if os.path.isfile(self.merge_intermdiate):
            os.remove(self.merge_intermdiate)

        perc_duplicates = None
        if os.path.isfile(self.picard_log):
            with open(self.picard_log) as log_fh:
                for line in log_fh:
                    if line.startswith("LIBRARY"):
                        fields = line.strip().split("\t")
                        duplication_idx = fields.index("PERCENT_DUPLICATION")
                        break
                fields = log_fh.next().strip().split("\t")
                perc_duplicates = float(fields[duplication_idx]) * 100
        elif os.path.isfile(self.dragen_log):
            with open(self.dragen_log) as d:
                for line in d:
                    line = line.strip()
                    if line.startswith("Number of duplicate reads"):
                        perc_duplicates = float(line.split("[")[-1].split("]")[0])
                        break
                    elif (line.startswith("MAPPING/ALIGNING SUMMARY") and
                          "Number of duplicate reads (marked)" in line):
                        perc_duplicates = float(line.split(" ")[-1])
                        break

        if perc_duplicates:
            with open(self.duplicates_file, "w") as out:
                out.write(str(perc_duplicates) + "\n")
        else:
            self.set_dsm_status(3) 
            raise ValueError("Could not find duplicate metrics in dragen log ({})".format(self.dragen_log))

        self.set_dsm_status(2) 

    # def requires(self):
        # return self.clone(CombineVariants), self.clone(PrintReads)
        # return self.clone(ValidateBAM)
        # return ValidateBAM(bam=self.scratch_bam, pseudo_prepid=self.pseudo_prepid, check_counts=self.sample_type != "CUSTOM_CAPTURE")#, FileExists(self.dragen_log)

class VariantCallingMetrics(GATKFPipelineTask):
    def pre_shell_commands(self):

        self.set_dsm_status(0) 

        self.metrics_raw = os.path.join(
            self.scratch_dir, self.sample_name + ".raw")
        self.metrics_raw_summary = self.metrics_raw + ".variant_calling_summary_metrics"
        self.metrics_raw_detail = self.metrics_raw + "variant_calling_detail_metrics"
        self.metrics_file = os.path.join(
            self.scratch_dir, self.name_prep + ".variant_calling_summary_metrics.txt")
        self.metrics_cmd = self.format_string(
            "{java} -XX:ParallelGCThreads=1 -jar {picard} CollectVariantCallingMetrics "
            "INPUT={annotated_vcf_gz} OUTPUT={metrics_raw} DBSNP={dbSNP} >& {log_file}")
        self.parser_cmd = self.format_string(
            "grep -v '^#' {metrics_raw_summary} | awk -f {transpose_awk} "
            "> {metrics_file}")
        self.shell_options.update(
            {"stdout":None, "stderr":None, "shell":True})
        self.commands = [self.metrics_cmd, self.parser_cmd]

    def requires(self):
        return self.clone(AnnotateVCF)

class ContaminationCheck(GATKFPipelineTask):
    """ Run VerifyBamID to check for sample contamination """
    def pre_shell_commands(self):

        self.set_dsm_status(0) 

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
    prept_completed_message = "GATK Pipeline" # Update the fact that the
    # pipeline has succeeded upon completion of this step

    def pre_shell_commands(self):

        self.set_dsm_status(0) 

        ########## dqm doesn't yet exist...

        ## Generic query to be used for updates 
        self.update_statement = """
        UPDATE {table} SET {field} = {value}
        WHERE pseudo_prepid = {pseudo_prepid}"""
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
        # self.master_table = "seqdbClone"
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

        # unbelieavable_wgs_hack=False

        if self.sample_type == "GENOME_AS_FAKE_EXOME":
            if not os.path.exists("/nfs/seqscratch_ssd/dsth/alignstats/PostReleaseMerge_PrePipelineEntry/{}.{}.txt".format(self.pseudo_prepid,self.sample_name)):
                hack_for_incorrect_single_rg_wgs="/nfs/goldstein/software/informatics/bin/alignstats -q 10 -i /nfs/seqscratch_ssd/ALIGNMENT/BUILD37/DRAGEN/GENOME_AS_FAKE_EXOME/{1}.{0}/{1}.{0}.bam \
                    -t /nfs/seqscratch_ssd/PIPELINE_DATA/ccds_regions.bed  -o /nfs/seqscratch_ssd/dsth/alignstats/PostReleaseMerge_PrePipelineEntry/{0}.{1}.txt".format(self.pseudo_prepid,self.sample_name)
                print("need to generate metrics for 'single rg' wgs - i.e. these don't really exist")
                if os.system(hack_for_incorrect_single_rg_wgs)!=0:
                #### single RG wgs samples - i.e. DON'T ACTUALLY EXISTS
                    raise ValueError("unable to generate metrics for single rg wgs sample")
                # unbelieavable_wgs_hack=True

        try:
            self.db = get_connection("seqdb")
            self.cur = self.db.cursor()

            self.cur.execute("select d.sample_name,s.chgvid,s.origid,d.is_external,p.externaldata,capture_kit from dragen_sample_metadata d "
              # + "join prepT p on d.pseudo_prepid=p.p_prepid "
              + "join Experiment e on d.experiment_id=e.id "
              + "join SampleT s on e.sample_id=s.sample_id "
              + "join prepT p on p.experiment_id=d.experiment_id "
              + "where d.pseudo_prepid = {}".format(self.pseudo_prepid))
            out = self.cur.fetchall()
            print(out)
            if ( len(out)!=1 and out[0][4] != None) or out[0][0]!=out[0][1]: 
                if out[0][5] == 'MCDMTOR':
                    print("i just don't care anymore")
                else:
                    raise ValueError("what is going on? : {}".format(out))
            vcf_sample_name = self.sample_name

            # dealing with messed up external samples as usual - the manifests are wrong!?!
            ##### hack for some broken external samples that were loaded incorrectly into LIMS
            if int(out[0][3])>1 and int(out[0][3])==4466:
                print("these are a real mess - we need to pull the key as is in the vcf but that gets truncated here anyway and it's already been purged from sampelt?!??!?")
                # from pprint import pprint as pp; pp(out)
                # pp(vars(self))
                ############# this is so odd. non-deterministic behaviour with cut?!?
                ############# just use vcf_format_name.split("\t")
                cmd='zcat {} | head -1000 | grep CHROM'.format(self.annotated_vcf_gz)
                # cmd='zcat {} | head -1000 | grep CHROM | cut -f10'.format(self.annotated_vcf_gz)
                wtf=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);
                vcf_format_name = wtf.stdout.readline().rstrip()
                vcf_format_name = vcf_format_name.split("\t")[9]
                print(cmd)
                print("from vcf we get '{}' : '{}'".format(vcf_format_name,self.annotated_vcf_gz))
                if vcf_format_name=="":
                    print("this is really creepy!?!")
                    time.sleep(3);
                    wtf=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);
                    vcf_format_name = wtf.stdout.readline().rstrip()
                    vcf_format_name = vcf_format_name.split("\t")[9]
                    if vcf_format_name=="":
                        print("wtf!?!")
                        os._exit(1) 
                    print("from vcf we get '{}' : '{}'".format(vcf_format_name,self.annotated_vcf_gz))
                    time.sleep(10);

                    ################ if this is still annoying simply extract the key directly since it's a terrible hack that turns off the check anyway!?!
                    ################ if this is still annoying simply extract the key directly since it's a terrible hack that turns off the check anyway!?!
                    ################ if this is still annoying simply extract the key directly since it's a terrible hack that turns off the check anyway!?!

                # vcf_format_name = subprocess.Popen('zcat {} | head -1000 | grep CHROM'.format(self.final_vcf),shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                ##### the names were fucked up in the manifest so cannot use that!?!
                ################################ all this effectively turns off any checks anyway so should just grab the keys and use it directly...
                if '.' in vcf_format_name:
                    vcf_sample_name=vcf_format_name.split('.')[0] # something removes suffix following .?!?
                elif '-' in vcf_format_name or '_' in vcf_format_name:
                    vcf_sample_name=vcf_format_name
                print("final name = {}".format(vcf_sample_name))
                    # vcf_sample_name=out[0][2].split('.')[0] # seems gatk or something removes suffix following .?!?
                    # print("updating key to '{}'".format(vcf_sample_name))
            ##### sn=BELATWEP278501 vs. here=BELATWEP2785 vs., 
            ##### > zcat /nfs/seqscratch_ssd/ALIGNMENT/BUILD37/DRAGEN/EXOME/BELATWEP278501.190927/BELATWEP278501.190927.analysisReady.vcf.gz | head -1000 | grep CHRO | cut -f10
            ##### BELATWEP2785.01
            ##### USANCH8771203 vs. USANCH8771-203

            # if int(out[0][3])>1 and int(out[0][3])!=3439 and out[0][0]!=out[0][2] and out[0][4] != None: # ignoring internal pushed through this way
                # vcf_sample_name=out[0][2].split('.')[0] # seems gatk or something removes suffix following .?!?

            with open(self.log_file, "w") as self.log_fh:
                self.get_variant_counts(vcf_sample_name)
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
                self.update_contamination_metrics()
                self.update_seqgender()
                self.update_qc_message()
                ## Update alignseqfile location here to keep track of scratch location
                # Due to bad decision the data location does not include the
                # sample's name and prep ID
                update_alignseqfile(
                    self.cur, os.path.dirname(self.scratch_dir), self.pseudo_prepid)
            self.db.commit()
            #for tmp_file in self.files_to_remove:
            #    for fn in glob(tmp_file):
            #        os.remove(fn)
        finally:
            if self.db.open:
                self.db.close()

        ######## I HATE ALL THIS NONSENSE SO MUCH THAT I JUST DON'T CARE ATM..
        if self.sample_type == "GENOME_AS_FAKE_EXOME":

            ###### no need for locking anything anymore!?!
            ########## just gonna hack it for now and integrate the newer pre-release and release intermediates in a completely non-systematic manner (these were new features patched on as separate pipelines...)
            gvcf_lazy=""
            metrics_lazy="/nfs/seqscratch_ssd/dsth/alignstats/PostReleaseMerge_PrePipelineEntry/{}.{}.txt".format(self.pseudo_prepid,self.sample_name)
            if self.sample_name[0:6] != "sqcudn" and self.sample_name!="3fcmt2":
                gvcf_lazy="/nfs/informatics/production/gvcf/{}.{}/{}.{}.gvcf.gz.md5sum".format(self.sample_name,self.pseudo_prepid,self.sample_name,self.pseudo_prepid)
                if not os.path.exists(gvcf_lazy):
                    raise ValueError("we're missing the full wgs gvcf md5 file ({})".format(gvcf_lazy))
                gvcf_lazy=gvcf_lazy[:-7]
                if not os.path.exists(gvcf_lazy):
                    raise ValueError("we're missing the full wgs gvcf file ({})".format(gvcf_lazy))
                print("\n\nusing '{}'".format(gvcf_lazy))
            print("\n\nusing '{}'".format(metrics_lazy))
            if not os.path.exists(metrics_lazy):
                #### single RG wgs samples - i.e. DON'T ACTUALLY EXISTS
                raise ValueError("we're missing the wgs metrics file ({})".format(metrics_lazy))
            r=re.compile("(WgsCoverageMe\w+)\": ([^,]+),")
            mean=None # mean=int()
            median=None # median=float()
            with open(metrics_lazy) as fh:
                for line in fh:
                    z=r.search(line) # if r.match(line): #### match is full string match!?!
                    if z: # print(z.groups())
                        if z.group(1)=="WgsCoverageMean":
                            mean=z.group(2)
                        elif z.group(1)=="WgsCoverageMedian":
                            median=z.group(2)
            print ("finally have mean={} and median={}".format(mean,median))
            q="update dragen_qc_metrics set WGS_Dragen_gVCF = '{}', WGS_Mean_Cov = {}, WGS_Median_Cov = {} where pseudo_prepid = {}".format(
            gvcf_lazy,mean,median,self.pseudo_prepid
            )
            print("final update will be '{}'".format(q))
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(q)
            print("updated = '{}' entries".format(cur.rowcount))
            if self.sample_name[0:6] != "sqcudn":
                if cur.rowcount!=1: # and unbelieavable_wgs_hack==False:
                    raise ValueError("unable to update wgs entries for sample ({})!".format(self.sample_name[0:6]))
                else:
                    print("update wgs metrics!")
            else:
                print("please get rid of this appling hack for sqc single rg wgs")
            db.commit()
            # print("\n\n\nTHIS IS DISABLED UNTIL INTEGRATE SOME OF THE ADDITIONAL OUT OF PIPE WGS METRICS LATER TODAY - this is now much simplified as no longer back patching - just grab the metrics but perhaps for lazyness just grab the relevant bits from the external scripts that were used?!?\n\n\n")
            # os._exit(1) 
            # raise ValueError("THIS IS DISABLED UNTIL INTEGRATE SOME OF THE ADDITIONAL OUT OF PIPE WGS METRICS LATER TODAY - this is now much simplified as no longer back patching - just grab the metrics but perhaps for lazyness just grab the relevant bits from the external scripts that were used?!?")

################ ProgrammingError: (1064, "You have an error in your SQL syntax; check the manual that corresponds to your MySQL server version for the right syntax to use near '?\n        WHERE pseudo_prepid = 42340' at line 1")
    def execute_query(self, query):
        self.log_fh.write(query + "\n")
        self.cur.execute(query)

    def add_sample_to_qc(self):
        """
        Initialize the sample in the qc table, use update statements later to update qc metrics
        """
        self.execute_query(self.format_string("""
        SELECT 1 FROM {qc_table} WHERE pseudo_prepid = {pseudo_prepid}"""))
        row = self.cur.fetchone()
        if not row:
            self.execute_query(self.format_string(
                "INSERT INTO {qc_table} (pseudo_prepid) VALUE ({pseudo_prepid})"))

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
        qc_failures = self.check_qc()
        if qc_failures:
            message.append(self.failed_qc)
            if self.issue_contamination_warning: ## Check for contamination warning
                message.append("Warning sample contamination is high, but below qc fail threshold")
            message.extend(qc_failures.values())
        else:
            message.append(self.passed_qc)
        final_message = '"{}"'.format(';'.join(message))
        self.update_database(self.qc_table,'QCMessage',final_message)

    def update_database(self, table, field, value):
        if value == '?':
            value = float("0.0")
        self.execute_query(self.format_string(
            self.update_statement, table=table, field=field, value=value))

    def get_metrics(self, query):
        try:
            self.cur.execute(query)
            db_val = self.cur.fetchall()
            return db_val
        except:
            raise ValueError(query)

    def update_alignment_metrics(self):
        self.alignment_metrics_map = {"TOTAL_READS":qc_metrics().total_reads,
                                      "PCT_PF_READS_ALIGNED":qc_metrics().pct_reads_aligned,
                                      "PF_MISMATCH_RATE":qc_metrics().pct_mismatch_rate}

        with open(self.alignment_metrics) as alignment_metrics_fh:
            first = True
            offset=0
            for line in alignment_metrics_fh:

                contents = line.strip().split(" ")

                if first:
                    first = False
                    if len(contents) == 2 and contents[-1] != "UNPAIRED":
                        raise ValueError("this should be simple unpaired reads")
                    elif len(contents) == 5 and ( contents[-2] != "PAIR" or contents[-1] != "UNPAIRED" ):
                        raise ValueError("this should be simple paired with unpaired reads (we should 'really' sum the paired&unpaired)")
                    elif len(contents) == 4 and contents[-1] != "PAIR":
                        raise ValueError("this should be simple paired only reads")

                if len(contents)==4 or len(contents)==5 or len(contents)==2:
                    field = contents[0]
                    value = contents[-2 if len(contents) == 5 else -1]

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

##### File "/nfs/goldstein/software/dragen_pipe/master/dragen/python/gatk_pipe.py", line 2004, in update_coverage_metrics
        with open(metrics_file) as metrics_fh:
            for line in metrics_fh:
                contents = line.strip().split(" ")
                if len(contents) > 1:
                    field, value  = contents[:2]
                    if field in metrics_hash:
                        db_field = metrics_hash[field]
                        if value == "NaN": # if math.isnan(value):
                            value="0.0"
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
################ IOError: [Errno 2] No such file or directory: '/nfs/seqscratch_ssd/ALIGNMENT/BUILD37/DRAGEN/CUSTOM_CAPTURE/nimhscz04C30271og3.15069/nimhscz04C30271og3.15069.variant_calling_summary_metrics.txt'
        with open(self.variant_call_out) as variant_call_metrics_fh:
            for line in variant_call_metrics_fh:
                field, value = line.strip().split(" ")[:2]
                if field in variant_call_parse:
                    db_field = variant_call_parse[field]
                    self.update_database(self.qc_table, db_field, value)
                if value == '?':
                    value = float("0.0")
                temp[field] = float(value)

        overall_titv = (
            (temp["NOVEL_SNPS"] * temp["NOVEL_TITV"] + temp["NUM_IN_DB_SNP"] * temp["DBSNP_TITV"]) /
            temp["TOTAL_SNPS"])
        self.update_database(self.qc_table, qc_metrics().titv, overall_titv)

        # homhet_ratio = (float(self.all_snv_hom + self.all_indel_hom) /
                        # (self.all_snv_het + self.all_indel_het))
        homhet_ratio=0.0
        if self.all_snv_het + self.all_indel_het > 0.0:
            homhet_ratio = (float(self.all_snv_hom + self.all_indel_hom) /
                            (self.all_snv_het + self.all_indel_het))

        self.update_database(self.qc_table, qc_metrics().homhet_ratio, homhet_ratio)

        snv_homhet_ratio = 0.0
        if self.all_snv_het > 0.0:
            snv_homhet_ratio = float(self.all_snv_hom) / self.all_snv_het
        self.update_database(self.qc_table, qc_metrics().snv_homhet_ratio, snv_homhet_ratio)

        if self.all_indel_het == 0.0:
            indel_homhet_ratio = 100.0
        else:
            indel_homhet_ratio = float(self.all_indel_hom) / self.all_indel_het
        self.update_database(self.qc_table, qc_metrics().indel_homhet_ratio, indel_homhet_ratio)

        x_homhet_ratio = 0.0
        if self.X_snv_het==0 and self.X_indel_het==0:
            if self.X_snv_hom>0 or self.X_indel_hom>0:
                x_homhet_ratio = 1.0
        else:
            x_homhet_ratio=(float(self.X_snv_hom + self.X_indel_hom) / (self.X_snv_het + self.X_indel_het))
                

        self.update_database(self.qc_table, qc_metrics().x_homhet_ratio, x_homhet_ratio)

    def update_contamination_metrics(self):
################ IOError: [Errno 2] No such file or directory: '/nfs/seqscratch_ssd/ALIGNMENT/BUILD37/DRAGEN/CUSTOM_CAPTURE/MZ200809564.5425/MZ200809564.5425.contamination.selfSM.txt'
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

        if self.sample_type in ("EXOME", "GENOME", "GENOME_AS_FAKE_EXOME"):
            query = """
            SELECT {cvg_type}
            FROM {qc_table}
            WHERE pseudo_prepid = {pseudo_prepid}"""
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
                    # https://redmine.igm.cumc.columbia.edu/issues/3022?issue_count=39&issue_position=24&next_issue_id=3179&prev_issue_id=3099
                    if isclose(0.0, Y_cvg) or ratio > 30:   # if isclose(0.0, Y_cvg) or ratio > 5:
                        seq_gender = "F"
                    elif ratio < 5:                         # elif ratio < 2:
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

            ratio = 1.0
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

        self.update_database(
            self.qc_table, 'SeqGender', '"{}"'.format(seq_gender))

    def check_alignment_rate(self):
        """
        Checks if the alignment metrics matches the appropriate thresholds
        Currently using only the perc_reads_aligned
        """

        query = self.format_string(
            """SELECT {pct_reads_aligned}
            FROM {qc_table}
            WHERE pseudo_prepid = {pseudo_prepid}""")
        result = self.get_metrics(query)
        from pprint import pprint as pp
        print("using dumbass query to retrieve transaction results?!? = {}".format(query))
        pp(result)
        perc_reads_aligned = float(result[0][0])
        return (perc_reads_aligned >= 0.70)

    def check_coverage(self):
        """
        Checks if the coverage metrics matches the appropriate thresholds
        """

        query = self.format_string(
            """SELECT {pct_bases5X}, {mean_cvg}
            FROM {qc_table}
            WHERE pseudo_prepid = {pseudo_prepid}""")
        result = self.get_metrics(query)
        pct_bases_5X = float(result[0][0])
        mean_cvg = float(result[0][1])
        return (pct_bases_5X >= 90 and mean_cvg >= 25)

    def check_duplicates(self):
        """
        Checks if the duplicate metrics matches the appropriate thresholds
        Currently using only the perc_duplicate_reads
        """
        query = self.format_string(
            """SELECT {pct_duplicate_reads}
            FROM {qc_table}
            WHERE pseudo_prepid = {pseudo_prepid}""")
        result = self.get_metrics(query)
        perc_duplicate_reads = float(result[0][0])
        if self.sample_type == "GENOME":
            return (perc_duplicate_reads <= 20)

        else:
            return (perc_duplicate_reads <= 30)

    def check_variant_calling(self):
        if self.sample_type == "GENOME":
            query = self.format_string(
                """SELECT {total_snps}
                FROM {qc_table} WHERE pseudo_prepid = {pseudo_prepid}""")
            result = self.get_metrics(query)
            if result[0][0] < 3000000:
                return False
        return not bool(check_vcf(self.vcf, self.check_variant_counts))
        # this "check" is really dumb and now deprecated
        #"""
        #Check if a minimal number of SNVs have been called, per sequencing type
        #Could very well fail on custom capture
        #"""
        #query = self.format_string(
        #    """SELECT {total_snps}
        #    FROM {qc_table} WHERE pseudo_prepid = {pseudo_prepid}""")
        #result = self.get_metrics(query)
        #num_snvs = int(result[0][0])

        #if self.sample_type == "GENOME":
        #    return num_snvs > 3000000
        #else:
        #    return num_snvs > 100000
        ## arguably should use something else for CUSTOM_CAPTURE...

    def check_contamination(self):
        """
        Checks if the contamination metrics matches the appropriate thresholds
        """
        query = self.format_string(
            """SELECT {contamination_value}
            FROM {qc_table}
            WHERE pseudo_prepid = {pseudo_prepid}""")
        result = self.get_metrics(query)
        contamination_value = float(result[0][0])
        if contamination_value < 0.08:
            if contamination_value >= 0.05:
                self.issue_contamination_warning = True
            return True
        elif contamination_value >= 0.08:
            return False

    def get_variant_counts(self,vcf_sample_name):
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
            header_fields[9] = header_fields[9].split(".")[0] 
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
                    variant_type = ("snv" if lREF == len(fields["ALT"]) else "indel")
                    # from pprint import pprint as pp; pp(fields)
                    # try:
                        # call_dict = dict(zip(fields["FORMAT"].split(":"),fields[self.sample_name].split(":")))
                    # except:
                        # raise ValueError("\n\nproblem with '{}' : '{}'\n".format(self.annotated_vcf_gz,fields,fields["FORMAT"],fields))
                    call_dict = dict(zip(fields["FORMAT"].split(":"),fields[vcf_sample_name].split(":")))
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
        yield self.clone(SplitAndSubsetDPBins)
        yield self.clone(AlignmentMetrics)
        yield self.clone(RunCvgMetrics)
        # yield self.clone(EntryChecks)
        yield self.clone(VariantCallingMetrics)
        yield self.clone(ContaminationCheck)

if __name__ == "__main__":
    os.environ['NO_PROXY'] = 'igm-atav01.igm.cumc.columbia.edu:8082'
    sys.exit(luigi.run())
