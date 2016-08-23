#!/nfs/goldstein/software/python2.7.7/bin/python2.7

import os
import shlex
import sys
import subprocess
import luigi
import MySQLdb
from luigi.contrib.sge import SGEJobTask
from dragen_globals import *
from dragen_sample import dragen_sample

"""
Run samples through a luigized GATK pipeline after finishing the
Dragen based alignment
"""

#Pipeline programs
gatk="/nfs/goldstein/software/GATK-3.6.0-nightly-2016-08-10-g9a77889/"
java="/nfs/goldstein/software/jdk1.8.0_05/bin/java"
bgzip="/nfs/goldstein/software/bin/bgzip"
tabix="/nfs/goldstein/software/bin/tabix"
picard="/nfs/goldstein/software/picard-tools-1.131/picard.jar"
bedtools="/nfs/goldstein/software/bedtools-2.25.0/bin/bedtools"
snpEff="/nfs/goldstein/software/snpEff/4.1/snpEff.jar"

#GATK parameters
max_mem="25"
#ref="/nfs/goldsteindata/refDB/HS_Build37/BWA_INDEX_hs37d5_BWAmem/hs37d5.fa"
ref="/scratch/HS_Build37/BWA_INDEX_hs37d5/hs37d5.fa"
hapmap="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/hapmap_3.3.b37.vcf"
omni="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/1000G_omni2.5.b37.vcf"
g1000="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/1000G_phase1.indels.b37.vcf"
Mills1000g="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/Mills_and_1000G_gold_standard.indels.b37.vcf"
dbSNP="/nfs/goldsteindata/refDB/dbSNP/dbsnp_147.b37.vcf"
interval="/nfs/goldstein/software/dragen_pipe/dragen/conf/hs37d5.intervals"

cfg = get_cfg()

class CopyBam(SGEJobTask):
    """class for copying dragen aligned bam to a scratch location"""
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"
    #job_name = self.sample_name + '.' + self.__class__.__name__ 

    def __init__(self, *args, **kwargs):
        super(CopyBam, self).__init__(*args, **kwargs)
        self.bam = "{base_directory}/{sample_name}/{sample_name}.bam".format(
            base_directory=self.base_directory,sample_name=self.sample_name)
        self.bam_index = "{base_directory}/{sample_name}/{sample_name}.bam.bai".format(
            base_directory=self.base_directory,sample_name=self.sample_name)
        self.scratch_bam = "{scratch}/{sample_name}/{sample_name}.bam".format(
            scratch=self.scratch,sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)
        db = get_connection("seqdb")
        try:
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
        db = get_connection("seqdb")
        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
        finally:
            if db.open:
                db.close()

        os.system("mkdir -p {0}/scripts".format(self.scratch + self.sample_name))
        os.system("mkdir -p {0}/logs".format(self.scratch + self.sample_name))
        cmd = ("rsync -a --timeout=20000 -r {bam} {bam_index} {scratch}/{sample_name}/"
              ).format(bam=self.bam,
                       scratch_bam=self.scratch_bam,
                       bam_index=self.bam_index,
                       scratch=self.scratch,
                       sample_name=self.sample_name)
        with open(self.script,'w') as o:
            o.write(cmd + "\n")
        subprocess.check_call(shlex.split(cmd))

        db = get_connection("seqcb")
        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
        finally:
            if db.open:
                db.close()

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)
        """
        yield luigi.LocalTarget("{scratch}/{sample_name}/{sample_name}.bam.bai".format(
            scratch=self.scratch,sample_name=self.sample_name))
        """
class RealignerTargetCreator(SGEJobTask):
    """class for creating targets for indel realignment BAMs from Dragen"""

    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')
    dbSNP = luigi.Parameter(default=dbSNP,
        description = 'dbSNP location')
    Mills1000g = luigi.Parameter(default=Mills1000g,
        description = 'Mills, Devin curated dataset')

    n_cpu = 4
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(RealignerTargetCreator, self).__init__(*args, **kwargs)
        self.interval_list = "{scratch}/{sample_name}/{sample_name}.interval_list".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)
        self.scratch_bam = "{scratch}/{sample_name}/{sample_name}.bam".format(
            scratch=self.scratch,sample_name=self.sample_name)

    def work(self):

        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T RealignerTargetCreator "
            "-L {interval} "
            "-I {scratch_bam} "
            "-o {interval_list} "
            "-known {Mills1000g} "
            "-known {dbSNP} "
            "-nt 4 ").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
                scratch_bam=self.scratch_bam,
                interval_list=self.interval_list,
                Mills1000g=Mills1000g,
                dbSNP=dbSNP)
        if not os.path.isdir(os.path.dirname(self.script)):
            os.makedirs(os.path.dirname(self.script))
        with open(self.script,'w') as o:
            o.write(cmd + "\n")
        subprocess.check_call(shlex.split(cmd))

    def requires(self):
        return self.clone(CopyBam)

    def output(self):
        yield luigi.LocalTarget(self.interval_list)


class IndelRealigner(SGEJobTask):
    """class to create BAM with realigned BAMs"""
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')
    dbSNP = luigi.Parameter(default=dbSNP,
        description = 'dbSNP location')
    Mills1000g = luigi.Parameter(default=Mills1000g,
        description = 'Mills, Devin curated dataset')

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(IndelRealigner, self).__init__(*args, **kwargs)
        self.interval_list = "{scratch}/{sample_name}/{sample_name}.interval_list".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.realn_bam = "{scratch}/{sample_name}/{sample_name}.realn.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)
        self.scratch_bam = "{scratch}/{sample_name}/{sample_name}.bam".format(
            scratch=self.scratch,sample_name=self.sample_name)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T IndelRealigner "
            "-L {interval} "
            "-I {scratch_bam} "
            "-o {realn_bam} "
            "-targetIntervals {interval_list} "
            "-maxReads 10000000 "
            "-maxInMemory 450000 "
            "-known {Mills1000g} "
            "-known {dbSNP}").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
                scratch_bam=self.scratch_bam,
                realn_bam=self.realn_bam,
                interval_list=self.interval_list,
                Mills1000g=Mills1000g,
                dbSNP=dbSNP)
        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))

    def requires(self):
        return self.clone(RealignerTargetCreator)

    def output(self):
        return luigi.LocalTarget(self.realn_bam)

class BaseRecalibrator(SGEJobTask):
    """class to create a recalibration table with realigned BAMs"""
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')
    dbSNP = luigi.Parameter(default=dbSNP,
        description = 'dbSNP location')
    Mills1000g = luigi.Parameter(default=Mills1000g,
        description = 'Mills, Devin curated dataset')
    n_cpu = 4
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(BaseRecalibrator, self).__init__(*args, **kwargs)
        self.bam = "{base_directory}/{sample_name}/{sample_name}.bam".format(
            base_directory=self.base_directory, sample_name=self.sample_name)
        self.realn_bam = "{scratch}/{sample_name}/{sample_name}.realn.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.recal_table = "{scratch}/{sample_name}/{sample_name}.recal_table".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)


    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T BaseRecalibrator "
            "-L {interval} "
            "-I {realn_bam} "
            "-nct 4 "
            "-o {recal_table} "
            "-L {capture_kit_bed} "
            "-knownSites {Mills1000g} "
            "-knownSites {dbSNP}").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
                realn_bam=self.realn_bam,
                recal_table=self.recal_table,
                Mills1000g=Mills1000g,
                capture_kit_bed=self.capture_kit_bed,
                dbSNP=dbSNP)
        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))


    def requires(self):
        return self.clone(IndelRealigner)

    def output(self):
        return luigi.LocalTarget(self.recal_table)

class PrintReads(SGEJobTask):
    """class to create BAM with realigned and recalculated BAMs"""
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')

    n_cpu = 4
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(PrintReads, self).__init__(*args, **kwargs)
        self.realn_bam = "{scratch}/{sample_name}/{sample_name}.realn.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.recal_table = "{scratch}/{sample_name}/{sample_name}.recal_table".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.recal_bam = "{scratch}/{sample_name}/{sample_name}.realn.recal.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
        # --disable_indel_quals are necessary to remove BI and BD tags in the bam file
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T PrintReads "
            "-L {interval} "
            "-I {realn_bam} "
            "--disable_indel_quals "
            "-BQSR {recal_table} "
            "-o {recal_bam} "
            "-nct 4").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
                realn_bam=self.realn_bam,
                recal_table=self.recal_table,
                recal_bam=self.recal_bam)

        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))

    def requires(self):
        return self.clone(BaseRecalibrator)

    def output(self):
        return luigi.LocalTarget(self.recal_bam)

class HaplotypeCaller(SGEJobTask):
    """class to create gVCFs"""
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    java = luigi.Parameter(default=java,description = 'java version used')
    gatk = luigi.Parameter(default=gatk,description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,description = 'reference genome location')
    dbSNP = luigi.Parameter(default=dbSNP,description = 'dbSNP location')
    n_cpu = 4
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(HaplotypeCaller, self).__init__(*args, **kwargs)
        self.realn_bam = "{scratch}/{sample_name}/{sample_name}.realn.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.recal_table = "{scratch}/{sample_name}/{sample_name}.recal_table".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.scratch_bam = "{scratch}/{sample_name}/{sample_name}.bam".format(
            scratch=self.scratch,sample_name=self.sample_name)
        self.recal_bam = "{scratch}/{sample_name}/{sample_name}.realn.recal.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.gvcf = "{scratch}/{sample_name}/{sample_name}.g.vcf.gz".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.gvcf_index = "{scratch}/{sample_name}/{sample_name}.g.vcf.gz.tbi".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T HaplotypeCaller "
            "-L {interval} "
            "-I {recal_bam} "
            "-o {gvcf} "
            "-stand_call_conf 20 "
            "-stand_emit_conf 20 "
            "--emitRefConfidence GVCF "
            "-GQB 5 -GQB 15 -GQB 20 -GQB 60 "
            "--variant_index_type LINEAR "
            "--variant_index_parameter 128000 "
            "--dbsnp {dbSNP} "
            "-nct 4").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
                recal_bam=self.recal_bam,
                gvcf=self.gvcf,
                dbSNP=dbSNP)

        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))

        rm_cmd = ['rm',self.realn_bam,self.scratch_bam]
        #print rm_cmd
        subprocess.call(rm_cmd)

    def requires(self):
        return self.clone(PrintReads)

    def output(self):
        return luigi.LocalTarget(self.gvcf_index)

class GenotypeGVCFs(SGEJobTask):
    """class to perfrom variant calling from gVCFs"""
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')
    dbSNP = luigi.Parameter(default=dbSNP,
        description = 'dbSNP location')

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(GenotypeGVCFs, self).__init__(*args, **kwargs)
        self.recal_table = "{scratch}/{sample_name}/{sample_name}.recal_table".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.gvcf = "{scratch}/{sample_name}/{sample_name}.g.vcf.gz".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.vcf = "{scratch}/{sample_name}/{sample_name}.raw.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T GenotypeGVCFs "
            "-L {interval} "
            "-o {vcf} "
            "-stand_call_conf 20 "
            "-stand_emit_conf 20 "
            "-V {gvcf}").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
                gvcf=self.gvcf,
                vcf=self.vcf)

        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))

    def requires(self):
        return self.clone(HaplotypeCaller)

    def output(self):
        return luigi.LocalTarget(self.vcf)

class SelectVariantsSNP(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(SelectVariantsSNP, self).__init__(*args, **kwargs)
        self.recal_table = "{scratch}/{sample_name}/{sample_name}.recal_table".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.vcf = "{scratch}/{sample_name}/{sample_name}.raw.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.snp_vcf = "{scratch}/{sample_name}/{sample_name}.snp.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T SelectVariants "
            "-L {interval} "
            "-V {vcf} "
            "-selectType SNP "
            "-o {snp_vcf}").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
                snp_vcf=self.snp_vcf,
                vcf=self.vcf)

        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))


    def requires(self):
        return self.clone(GenotypeGVCFs)

    def output(self):
        return luigi.LocalTarget(self.snp_vcf)

class SelectVariantsINDEL(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(SelectVariantsINDEL, self).__init__(*args, **kwargs)
        self.recal_table = "{scratch}/{sample_name}/{sample_name}.recal_table".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.vcf = "{scratch}/{sample_name}/{sample_name}.raw.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_vcf = "{scratch}/{sample_name}/{sample_name}.indel.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T SelectVariants "
            "-L {interval} "
            "-V {vcf} "
            "-selectType INDEL "
            "-o {indel_vcf}").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
                indel_vcf=self.indel_vcf,
                vcf=self.vcf)

        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))


    def requires(self):
        if self.sample_type == 'exome':
            return self.clone(ApplyRecalibrationSNP)
        if self.sample_type == 'genome':
            return self.clone(ApplyRecalibrationSNP)
        elif self.sample_type == 'custom_capture':
            return self.clone(VariantFiltrationSNP)
        else:
            raise Exception, "Sample type: %s not supported in this module" % self.sample_type

    def output(self):
        return luigi.LocalTarget(self.indel_vcf)


class VariantRecalibratorSNP(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')
    dbSNP = luigi.Parameter(default=dbSNP,
        description = 'dbSNP location')

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(VariantRecalibratorSNP, self).__init__(*args, **kwargs)
        self.snp_vcf = "{scratch}/{sample_name}/{sample_name}.snp.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.snp_recal = "{scratch}/{sample_name}/{sample_name}.snp.recal".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.snp_tranches = "{scratch}/{sample_name}/{sample_name}.snp.tranches".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.snp_rscript = "{scratch}/{sample_name}/{sample_name}.snp.rscript".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T VariantRecalibrator "
            "-L {interval} "
            "--input {snp_vcf} "
            "-an DP "
            "-an QD "
            "-an FS "
            "-an SOR "
            "-an MQ "
            "-an MQRankSum "
            "-an ReadPosRankSum "
            "-mode SNP "
            "--maxGaussians 4 "
            "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 "
            "-recalFile {snp_recal} "
            "-tranchesFile {snp_tranches} "
            "-rscriptFile {snp_rscript} "
            "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} "
            "-resource:omni,known=false,training=true,truth=false,prior=12.0 {omni} "
            "-resource:1000G,known=false,training=true,truth=false,prior=10.0 {g1000} "
            "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbSNP} "
            ).format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
                snp_vcf=self.snp_vcf,
                snp_recal=self.snp_recal,
                snp_tranches=self.snp_tranches,
                snp_rscript=self.snp_rscript,
                dbSNP=dbSNP,
                hapmap=hapmap,
                omni=omni,
                g1000=g1000)
        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))

    def requires(self):
        return self.clone(SelectVariantsSNP)


    def output(self):
        #return luigi.LocalTarget(self.snp_rscript),luigi.LocalTarget(self.snp_tranches),luigi.LocalTarget(self.snp_recal)
        return luigi.LocalTarget(self.snp_recal)

        ####double check dbsnp version ####

class VariantRecalibratorINDEL(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')
    dbSNP = luigi.Parameter(default=dbSNP,
        description = 'dbSNP location')
    Mills1000g = luigi.Parameter(default=Mills1000g,
        description = 'Mills, Devin curated dataset')

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(VariantRecalibratorINDEL, self).__init__(*args, **kwargs)
        self.indel_vcf = "{scratch}/{sample_name}/{sample_name}.indel.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_recal = "{scratch}/{sample_name}/{sample_name}.indel.recal".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_rscript = "{scratch}/{sample_name}/{sample_name}.indel.rscript".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_tranches = "{scratch}/{sample_name}/{sample_name}.indel.tranches".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T VariantRecalibrator "
            "-L {interval} "
            "--input {indel_vcf} "
            "-an QD "
            "-an DP "
            "-an FS "
            "-an SOR "
            "-an MQRankSum "
            "-an ReadPosRankSum "
            "-mode INDEL "
            "--maxGaussians 4 "
            "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 "
            "-recalFile {indel_recal} "
            "-tranchesFile {indel_tranches} "
            "-resource:mills,known=true,training=true,truth=true,prior=12.0 {Mills1000g} "
            ).format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
                indel_vcf=self.indel_vcf,
                indel_recal=self.indel_recal,
                indel_tranches=self.indel_tranches,
                indel_rscript=self.indel_rscript,
                Mills1000g=Mills1000g)
        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))

    def requires(self):
        return self.clone(SelectVariantsINDEL)

    def output(self):
        return luigi.LocalTarget(self.indel_recal)

class ApplyRecalibrationSNP(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(ApplyRecalibrationSNP, self).__init__(*args, **kwargs)
        self.vcf = "{scratch}/{sample_name}/{sample_name}.snp.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.snp_recal = "{scratch}/{sample_name}/{sample_name}.snp.recal".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.snp_tranches = "{scratch}/{sample_name}/{sample_name}.snp.tranches".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.snp_filtered = "{scratch}/{sample_name}/{sample_name}.snp.filtered.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T ApplyRecalibration "
            "-L {interval} "
            "-input {vcf} "
            "-tranchesFile {snp_tranches} "
            "-recalFile {snp_recal} "
            "-o {snp_filtered} "
            "--ts_filter_level 99.0 "
            "-mode SNP ").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
                vcf=self.vcf,
                snp_tranches=self.snp_tranches,
                snp_filtered=self.snp_filtered,
                snp_recal=self.snp_recal)
        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))

    def requires(self):
        return self.clone(VariantRecalibratorSNP)

    def output(self):
        return luigi.LocalTarget(self.snp_filtered)

class ApplyRecalibrationINDEL(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(ApplyRecalibrationINDEL, self).__init__(*args, **kwargs)
        self.indel_vcf = "{scratch}/{sample_name}/{sample_name}.indel.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_recal = "{scratch}/{sample_name}/{sample_name}.indel.recal".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_tranches = "{scratch}/{sample_name}/{sample_name}.indel.tranches".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_filtered = "{scratch}/{sample_name}/{sample_name}.indel.filtered.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T ApplyRecalibration "
            "-L {interval} "
            "-input {indel_vcf} "
            "-tranchesFile {indel_tranches} "
            "-recalFile {indel_recal} "
            "-o {indel_filtered} "
            "--ts_filter_level 99.0 "
            "-mode INDEL ").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
                indel_vcf=self.indel_vcf,
                indel_recal=self.indel_recal,
                indel_tranches=self.indel_tranches,
                indel_filtered=self.indel_filtered)
       with open(self.script,'w') as o:
           o.write(cmd + "\n")
           subprocess.check_call(shlex.split(cmd))

    def requires(self):
        return self.clone(VariantRecalibratorINDEL)

    def output(self):
        return luigi.LocalTarget(self.indel_filtered)

class VariantFiltrationSNP(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(VariantFiltrationSNP, self).__init__(*args, **kwargs)
        self.snp_vcf = "{scratch}/{sample_name}/{sample_name}.snp.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.snp_filtered = "{scratch}/{sample_name}/{sample_name}.snp.filtered.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T VariantFiltration "
            "-L {interval} "
            "-V {snp_vcf} "
            "--filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" "
            "--filterName \"SNP_filter\" "
            "-o {snp_filtered} ").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
                snp_vcf=self.snp_vcf,
                snp_filtered=self.snp_filtered)
       with open(self.script,'w') as o:
           o.write(cmd + "\n")
           subprocess.check_call(shlex.split(cmd))

    def requires(self):
        return self.clone(SelectVariantsSNP)

    def output(self):
        return luigi.LocalTarget(self.snp_filtered)

class VariantFiltrationINDEL(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(VariantFiltrationINDEL, self).__init__(*args, **kwargs)
        self.indel_vcf = "{scratch}/{sample_name}/{sample_name}.indel.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_filtered = "{scratch}/{sample_name}/{sample_name}.indel.filtered.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T VariantFiltration "
            "-L {interval} "
            "-V {indel_vcf} "
            "--filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" "
            "--filterName \"INDEL_filter\" "
            "-o {indel_filtered}").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
                indel_vcf=self.indel_vcf,
                indel_filtered=self.indel_filtered)
       with open(self.script,'w') as o:
           o.write(cmd + "\n")
           subprocess.check_call(shlex.split(cmd))

    def requires(self):
        return self.clone(SelectVariantsINDEL)

    def output(self):
        return luigi.LocalTarget(self.indel_filtered)


class CombineVariants(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)

    java = luigi.Parameter(default=java,
        description = 'java version used')
    tabix = luigi.Parameter(default=tabix,
        description = 'tabix version used')
    picard = luigi.Parameter(default=picard,
        description = 'picard version used')
    bgzip = luigi.Parameter(default=bgzip,
        description = 'bgzip version used')

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(CombineVariants, self).__init__(*args, **kwargs)
        self.snp_filtered = "{scratch}/{sample_name}/{sample_name}.snp.filtered.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.vcf = "{scratch}/{sample_name}/{sample_name}.raw.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_filtered = "{scratch}/{sample_name}/{sample_name}.indel.filtered.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.final_vcf = "{scratch}/{sample_name}/{sample_name}.analysisReady.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.tmp_vcf = "{scratch}/{sample_name}/{sample_name}.tmp.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
        """Merges SNP and INDEL vcfs.  Using the filtered SNP vcf header as a
        base, variant type/sample type specific ##FILTERs are added to the header.
        Additionally the AnnoDBID annotation is added.  After finishing reading
        through the header the SNP vcf is read again with only variants
        being outputted.  The same happens with the INDEL vcf.  The resulting
        vcf is then sorted producing the analysisReady.vcf.
        """
        filter_flag = 0
        with open(self.tmp_vcf,'w') as vcf_out, open(self.snp_filtered) as header:
            for line in header.readlines():
                if line[0] == '#':
                    if line[0:8] == '##FILTER' and filter_flag == 0 :
                        filter_flag = 1
                        #SNP specific FILTER whether using VQSR or snp filtering
                        if self.sample_type == 'exome' or self.sample_type =='genome':
                            vcf_out.write('##FILTER=<ID=VQSRTrancheSNP99.00to99.90,Description="Truth sensitivity tranche level for SNP model\n')
                            vcf_out.write('##FILTER=<ID=VQSRTrancheSNP99.90to100.00+,Description="Truth sensitivity tranche level for SNP model\n')
                            vcf_out.write('##FILTER=<ID=VQSRTrancheSNP99.90to100.00,Description="Truth sensitivity tranche level for SNP model\n')
                        else:
                            vcf_out.write('##FILTER=<ID=SNP_filter,Description="QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0">\n')

                        #Indel specific filters whether using VQSR or indel filtering
                        if self.sample_type =='genome':
                            vcf_out.write('##FILTER=<ID=VQSRTrancheINDEL99.00to99.90,Description="Truth sensitivity tranche level for INDEL model\n')
                            vcf_out.write('##FILTER=<ID=VQSRTrancheINDEL99.90to100.00+,Description="Truth sensitivity tranche level for INDEL model\n')
                            vcf_out.write('##FILTER=<ID=VQSRTrancheINDEL99.90to100.00,Description="Truth sensitivity tranche level for INDEL model\n')
                        else:
                            vcf_out.write('##FILTER=<ID=INDEL_filter,Description="QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0">\n')

                    #AnnoDBID annotation will be added during the annotation pipeline
                    if line[0:13] == "##INFO=<ID=AN":
                            vcf_out.write('##INFO=<ID=AnnoDBID,Number=1,Type=String,Description="AnnoDBID">\n')
                    vcf_out.write(line)
                else:
                    break

            with open(self.snp_filtered) as snps:

                for snp in snps.readlines():
                    if snp[0] != '#':
                        vcf_out.write(snp)
            with open(self.indel_filtered) as indels:
                for indel in indels.readlines():
                    if indel[0] != '#':
                        vcf_out.write(indel)

        sort_cmd = ("{java} -jar {picard} "
            "SortVcf "
            "I={tmp_vcf} "
            "O={final_vcf}").format(java=java,
                picard=picard,
                tmp_vcf=self.tmp_vcf,
                final_vcf=self.final_vcf)

        with open(self.script,'w') as o:
            o.write(sort_cmd + "\n")

        subprocess.check_call(shlex.split(sort_cmd))

    def requires(self):
        if self.sample_type == 'exome':
            return self.clone(VariantFiltrationINDEL)
        elif self.sample_type == 'custom_capture':
            return self.clone(VariantFiltrationINDEL)
        elif self.sample_type == 'genome':
            return self.clone(ApplyRecalibrationINDEL)
        else:
            raise Exception, "Sample type: %s not supported in this module" % self.sample_type
    def output(self):
        return luigi.LocalTarget(self.final_vcf)

class AnnotateVCF(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    scratch = luigi.Parameter()
    sample_type = luigi.Parameter()

    java = luigi.Parameter(default=java,
        description = 'java version used')
    tabix = luigi.Parameter(default=tabix,
        description = 'tabix version used')
    bgzip = luigi.Parameter(default=bgzip,
        description = 'bgzip version used')
    snpEff = luigi.Parameter(default=snpEff,
        description = 'snpEff version used')

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(AnnotateVCF, self).__init__(*args, **kwargs)
        self.snp_filtered = "{scratch}/{sample_name}/{sample_name}.snp.filtered.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.final_vcf = "{scratch}/{sample_name}/{sample_name}.analysisReady.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.annotated_vcf = "{scratch}/{sample_name}/{sample_name}.analysisReady.annotated.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.annotated_vcf_gz = "{scratch}/{sample_name}/{sample_name}.analysisReady.annotated.vcf.gz".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.annotated_vcf_gz_index = "{scratch}/{sample_name}/{sample_name}.analysisReady.annotated.vcf.gz.tbi".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
        #SnpEff parameters
        snpeff_cfg="/nfs/goldstein/software/snpEff/4.1/snpEff.config"
        snpeff_options=""
        intervals="/nfs/goldstein/goldsteinlab/software/snpEff_3_3/data/GRCh37.73/intron_exon_boundaries.unique.custom_interval.bed"

        snpEff_cmd = ("{java} -Xmx5G -jar {snpEff} eff "
                     "GRCh37.74 -c {snpeff_cfg} "
                     "-interval {intervals} "
                     "-v -noMotif -noNextProt -noLog -nodownload -noStats "
                     "-o vcf {final_vcf}"
                     ).format(java=java,
                              snpEff=snpEff,
                              snpeff_cfg=snpeff_cfg,
                              intervals=intervals,
                              final_vcf=self.final_vcf)
        bgzip_cmd = ("{bgzip} {annotated_vcf}").format(bgzip=bgzip,annotated_vcf=self.annotated_vcf)
        tabix_cmd = ("{tabix} {annotated_vcf_gz}").format(tabix=tabix,annotated_vcf_gz=self.annotated_vcf_gz)

        with open(self.script,'w') as o:
            o.write(snpEff_cmd + "\n")
            o.write(bgzip_cmd + "\n")
            o.write(tabix_cmd + "\n")

        with open(self.annotated_vcf,'w') as vcf_out, \
            open(self.annotated_vcf + ".log", "w") as log_fh:
                p = subprocess.Popen(shlex.split(snpEff_cmd), stdout=vcf_out, stderr=log_fh)
                p.wait()
        if p.returncode:
            raise subprocess.CalledProcessError(p.returncode, cmd)

        subprocess.check_call(shlex.split(snpEff_cmd))
        #subprocess.check_call(shlex.split(bgzip_cmd))
        #subprocess.check_call(shlex.split(tabix_cmd))

    def requires(self):
        return self.clone(CombineVariants)

    def output(self):
        return luigi.LocalTarget(self.annotated_vcf_gz_index)

class ArchiveSample(SGEJobTask):
    """Archive samples on Amplidata"""

    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)
    #dragen_id = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    """
    db = MySQLdb.connect(db="sequenceDB", read_default_group="clientsequencedb",
            read_default_file="/nfs/goldstein/software/dragen/.my.cnf")
    curs = db.cursor()
    """

    def __init__(self, *args, **kwargs):
        super(ArchiveSample, self).__init__(*args, **kwargs)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)
        self.recal_bam = "{scratch}/{sample_name}/{sample_name}.realn.recal.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.recal_bam_index = "{scratch}/{sample_name}/{sample_name}.realn.recal.bai".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.bam = "{base_directory}/{sample_name}/{sample_name}.bam".format(
            base_directory=self.base_directory, sample_name=self.sample_name)
        self.annotated_vcf_gz = "{scratch}/{sample_name}/{sample_name}.analysisReady.annotated.vcf.gz".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.g_vcf_gz = "{scratch}/{sample_name}/{sample_name}.g.vcf.gz".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.g_vcf_gz_index = "{scratch}/{sample_name}/{sample_name}.g.vcf.gz.tbi".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.annotated_vcf_gz_index = "{scratch}/{sample_name}/{sample_name}.analysisReady.annotated.vcf.gz.tbi".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.copy_complete = "{base_directory}/{sample_name}/copy_complete".format(
            base_directory=self.base_directory, sample_name=self.sample_name)
        self.script_dir = "{scratch}/{sample_name}/scripts".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.scratch_dir = "{scratch}/{sample_name}".format(
            scratch=self.scratch,sample_name=self.sample_name)

    def work(self):
        cmd = ("rsync -a --timeout=25000 -r "
              "{script_dir} {recal_bam} {recal_bam_index} {annotated_vcf_gz} "
              "{annotated_vcf_gz_index} {g_vcf_gz} {base_directory}/{sample_name} "
              "{g_vcf_gz_index}"
              ).format(recal_bam=self.recal_bam,
                      recal_bam_index=self.recal_bam_index,
                      script_dir=self.script_dir,
                      annotated_vcf_gz=self.annotated_vcf_gz,
                      annotated_vcf_gz_index=self.annotated_vcf_gz_index,
                      base_directory=self.base_directory,
                      g_vcf_gz=self.g_vcf_gz,
                      g_vcf_gz_index=self.g_vcf_gz_index,
                      sample_name=self.sample_name)
        with open(self.script,'w') as o:
            o.write(cmd + "\n")
        subprocess.check_call(shlex.split(cmd))
        subprocess.call(['touch',self.copy_complete])

        # Original dragen BAM could be technically deleted earlier after the
        # realigned BAM has been created on scratch space but it is safer to 
        # delete after the final realigned, recalculated BAM has been archived
        # since our scratch space has failed in the past.

        rm_cmd = ['rm',self.bam]
        rm_folder_cmd = ['rm','-rf',self.scratch_dir]
        #subprocess.call(rm_cmd)
        #subprocess.call(rm_folder_cmd)

    def requires(self):
        return self.clone(AnnotateVCF)

    def output(self):
        return luigi.LocalTarget(self.copy_complete)


class SQLTarget(luigi.Target):
    """Target describing verification of the entries in the database
    """
    def __init__(self, pseudo_prepid, pipeline_step_id):
        self.pseudo_prepid = pseudo_prepid
        self.pipeline_step_id = pipeline_step_id

    def exists(self):
        db = get_connection("seqdb")
        try:
            cur = db.cursor()
            cur.execute(GET_STEP_STATUS.format(
                pseudo_prepid=self.pseudo_prepid,
                pipeline_step_id=self.pipeline_step_id))
            row = cur.fetchone()
            if row:
                if row[0] == 1:
                    return True
                else:
                    return False
            else:
                cur.execute(INSERT_PIPELINE_STEP.format(
                    pseudo_prepid=self.pseudo_prepid,
                    pipeline_step_id=self.pipeline_step_id))
                db.commit()
                return False
        finally:
            if db.open:
                db.close()

