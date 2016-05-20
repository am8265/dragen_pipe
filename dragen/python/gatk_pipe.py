import luigi
import os
import shlex
import sys
import subprocess
from dragen_sample import dragen_sample
#from luigi import configuration
from luigi.contrib.sge import SGEJobTask, LocalSGEJobTask

"""
Run samples through a luigized GATK pipeline after they have finished the
Dragen based alignment
"""

gatk="/nfs/goldstein/software/GATK-3.5.0"
java="/nfs/goldstein/goldsteinlab/software/java/jdk1.7.0_03/bin/java"
max_mem="15"
ref="/nfs/goldsteindata/refDB/HS_Build37/BWA_INDEX_hs37d5_BWAmem/hs37d5.fa"
hapmap="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/hapmap_3.3.b37.vcf"
omni="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/1000G_omni2.5.b37.vcf"
g1000="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/1000G_phase1.indels.b37.vcf"
Mills1000g="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/Mills_and_1000G_gold_standard.indels.b37.vcf"
dbSNP="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/dbsnp_138.b37.vcf"


class RootTask(LocalSGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()

    def requires(self):
        if sample_type == 'genome':
            return self.clone(IndelRealigner)
    def work(self):
        os.system('touch /home/jb3816/{0}.complete'.format(self.sample_name))
    def output(self):
        return luigi.LocalTarget('/home/jb3816/{0}.complete'.format(self.sample_name))

class RealignerTargetCreator(LocalSGEJobTask):
    """class for creating targets for indel realignment BAMs from Dragen"""

    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
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

    def __init__(self, *args, **kwargs):
        super(RealignerTargetCreator, self).__init__(*args, **kwargs)
        self.bam = "{base_directory}/{sample_name}/{sample_name}.bam".format(
            base_directory=self.base_directory, sample_name=self.sample_name)
        self.interval_list = "{scratch}/{sample_name}/{sample_name}.interval_list".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.RealignerTargetCreator_script = "{scratch}/{sample_name}/scripts/RealignerTargetCreator.sh".format(
            scratch=self.scratch, sample_name=self.sample_name)

    def work(self):
        os.system("mkdir -p {0}/scripts".format(self.scratch + self.sample_name))
        os.system("mkdir -p {0}/logs".format(self.scratch + self.sample_name))

        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T RealignerTargetCreator "
            "-I {bam} "
            "-o {interval_list} "
            "-known {Mills1000g} "
            "-known {dbSNP} "
            "-nt 4 ").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                bam=self.bam,
                interval_list=self.interval_list,
                Mills1000g=Mills1000g,
                dbSNP=dbSNP)
        #os.system("echo {cmd}".format(cmd=cmd))
        with open(self.RealignerTargetCreator_script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))
        with self.output().open("w") as out:
            out.write("test")
        #os.system("touch {0}".format(self.interval_list))

    def output(self):
        return luigi.LocalTarget(self.interval_list)

class IndelRealigner(LocalSGEJobTask):
    """class to create BAM with realigned BAMs"""
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
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

    def __init__(self, *args, **kwargs):
        super(IndelRealigner, self).__init__(*args, **kwargs)
        self.bam = "{base_directory}/{sample_name}/{sample_name}.bam".format(
            base_directory=self.base_directory, sample_name=self.sample_name)
        self.interval_list = "{scratch}/{sample_name}/{sample_name}.interval_list".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.realn_bam = "{scratch}/{sample_name}/{sample_name}.realn.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T IndelRealigner "
            "-I {bam} "
            "-o {realn_bam} "
            "-targetIntervals {interval_list} "
            "-maxReads 10000000 "
            "-maxInMemory 450000 "
            "-known {Mills1000g} "
            "-known {dbSNP}").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                bam=self.bam,
                realn_bam=self.realn_bam,
                interval_list=self.interval_list,
                Mills1000g=Mills1000g,
                dbSNP=dbSNP)
        os.system("echo {0}".format(cmd))
        os.system("touch {0}".format(self.realn_bam))

    def requires(self):
        return self.clone(RealignerTargetCreator)

    def output(self):
        return luigi.LocalTarget(self.realn_bam)

class BaseRecalibrator(LocalSGEJobTask):
    """class to create a recalibration table with realigned BAMs"""
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
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

    def __init__(self, *args, **kwargs):
        super(BaseRecalibrator, self).__init__(*args, **kwargs)
        self.bam = "{base_directory}/{sample_name}/{sample_name}.bam".format(
            base_directory=self.base_directory, sample_name=self.sample_name)
        self.realn_bam = "{scratch}/{sample_name}/{sample_name}.realn.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.recal_table = "{scratch}/{sample_name}/{sample_name}.recal_table".format(
            scratch=self.scratch, sample_name=self.sample_name)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T BaseRecalibrator "
            "-I {realn_bam} "
            "-nct 4 "
            "-o {recal_table} "
            "-known {Mills1000g} "
            "-known {dbSNP}").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                realn_bam=self.realn_bam,
                recal_table=self.recal_table,
                Mills1000g=Mills1000g,
                dbSNP=dbSNP)

        os.system("echo {0}".format(cmd))
        os.system("touch {0}".format(self.recal_table))

        rm_cmd = ['rm',self.bam]
        print rm_cmd
        #subprocess.call(rm_cmd)

    def requires(self):
        return self.clone(IndelRealigner)

    def output(self):
        return luigi.LocalTarget(self.recal_table)

class PrintReads(LocalSGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')

    def __init__(self, *args, **kwargs):
        super(PrintReads, self).__init__(*args, **kwargs)
        self.realn_bam = "{scratch}/{sample_name}/{sample_name}.realn.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.recal_table = "{scratch}/{sample_name}/{sample_name}.recal_table".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.recal_bam = "{scratch}/{sample_name}/{sample_name}.realn.recal.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T PrintReads "
            "-I {realn_bam} "
            "-BQSR {recal_table} "
            "-o {recal_bam} "
            "-nct 4").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                realn_bam=self.realn_bam,
                recal_table=self.recal_table,
                recal_bam=self.recal_bam)

        os.system("echo {0}".format(cmd))
        os.system("touch {0}".format(self.recal_bam))

    def requires(self):
        return self.clone(BaseRecalibrator)

    def output(self):
        return luigi.LocalTarget(self.recal_bam)

class HaplotypeCaller(LocalSGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
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

    def __init__(self, *args, **kwargs):
        super(HaplotypeCaller, self).__init__(*args, **kwargs)
        self.recal_table = "{scratch}/{sample_name}/{sample_name}.recal_table".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.recal_bam = "{scratch}/{sample_name}/{sample_name}.realn.recal.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.vcf = "{scratch}/{sample_name}/{sample_name}.raw.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)

    def work(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T HaplotypeCaller "
            "-I {recal_bam} "
            "-o {vcf} "
            "-stand_call_conf 20 "
            "-stand_emit_conf 20 "
            "--dbsnp {dbSNP} "
            "-nct 4").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                recal_bam=self.recal_bam,
                vcf=self.vcf,
                dbSNP=dbSNP)

       rm_cmd = ['rm',self.recal_bam]
       print rm_cmd
       #subprocess.call(rm_cmd)
       os.system("echo {0}".format(cmd))
       os.system("touch {0}".format(self.vcf))

    def requires(self):
        return self.clone(PrintReads)

    def output(self):
        return luigi.LocalTarget(self.vcf)

class VariantRecalibratorSNP(LocalSGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
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

    def __init__(self, *args, **kwargs):
        super(VariantRecalibratorSNP, self).__init__(*args, **kwargs)
        self.vcf = "{scratch}/{sample_name}/{sample_name}.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.snp_recal = "{scratch}/{sample_name}/{sample_name}.snp.recal".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.snp_tranches = "{scratch}/{sample_name}/{sample_name}.snp.tranches".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.snp_rscript = "{scratch}/{sample_name}/{sample_name}.snp.rscript".format(
            scratch=self.scratch, sample_name=self.sample_name)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T VariantRecalibrator "
            "-I {vcf} "
            "-an DP "
            "-an QD "
            "-an FS "
            "-an SOR "
            "-an MQRankSum "
            "-an ReadPosRankSum "
            "-an MQ "
            "-mode SNP "
            "-tranche 100.0 "
            "-tranche 99.9 "
            "-tranche 99.0 "
            "-tranche 90.0 "
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
                vcf=self.vcf,
                snp_recal=self.snp_recal,
                snp_tranches=self.snp_tranches,
                snp_rscript=self.snp_rscript,
                dbSNP=dbSNP,
                hapmap=hapmap,
                omni=omni,
                g1000=g1000)
        os.system("echo {0}".format(cmd))
        os.system("touch {0}".format(self.snp_tranches))
        os.system("touch {0}".format(self.snp_rscript))
        os.system("touch {0}".format(self.snp_recal))

    def requires(self):
        return self.clone(HaplotypeCaller)

    def output(self):
        #return luigi.LocalTarget(self.snp_rscript),luigi.LocalTarget(self.snp_tranches),luigi.LocalTarget(self.snp_recal)
        return luigi.LocalTarget(self.snp_recal)

        ####double check resources version ####
        ####double check dbsnp version ####

class VariantRecalibratorINDEL(LocalSGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
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

    def __init__(self, *args, **kwargs):
        super(VariantRecalibratorINDEL, self).__init__(*args, **kwargs)
        self.vcf = "{scratch}/{sample_name}/{sample_name}.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_recal = "{scratch}/{sample_name}/{sample_name}.indel.recal".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_tranches = "{scratch}/{sample_name}/{sample_name}.indel.tranches".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_rscript = "{scratch}/{sample_name}/{sample_name}.indel.rscript".format(
            scratch=self.scratch, sample_name=self.sample_name)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T VariantRecalibrator "
            "-I {vcf} "
            "-tranche 100.0 "
            "-tranche 99.9 "
            "-tranche 99.0 "
            "-tranche 90.0 "
            "-recalFile {indel_recal} "
            "-tranchesFile {indel_tranches} "
            "-rscriptFile {indel_rscript} "
            "--maxGaussians 4 "
            "-resource:mills,known=false,training=true,truth=true,prior=12.0 {Mills1000g} "
            "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbSNP} "
            "-an QD "
            "-an DP "
            "-an FS "
            "-an SOR "
            "-an ReadPosRankSum "
            "-an MQRankSum "
            "-an InbreedingCoeff "
            "-mode INDEL ").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                vcf=self.vcf,
                indel_recal=self.indel_recal,
                indel_tranches=self.indel_tranches,
                indel_rscript=self.indel_rscript,
                dbSNP=dbSNP,
                Mills1000g=Mills1000g)
        os.system("echo {0}".format(cmd))
        os.system("touch {0}".format(self.indel_recal))
        os.system("touch {0}".format(self.indel_tranches))
        os.system("touch {0}".format(self.indel_rscript))

    def requires(self):
        return self.clone(VariantRecalibratorSNP)

    def output(self):
        return luigi.LocalTarget(self.indel_recal)

class ApplyRecalibrationSNP(LocalSGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')

    def __init__(self, *args, **kwargs):
        super(ApplyRecalibrationSNP, self).__init__(*args, **kwargs)
        self.vcf = "{scratch}/{sample_name}/{sample_name}.snp.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.snp_recal = "{scratch}/{sample_name}/{sample_name}.snp.recal2".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.snp_tranches = "{scratch}/{sample_name}/{sample_name}.snp.tranches".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.snp_filtered = "{scratch}/{sample_name}/{sample_name}.snp.filtered.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)

    def work(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T ApplyRecalibration "
            "-I {vcf} "
            "-tranchesFile {snp_tranches} "
            "-recalFile {snp_recal} "
            "-o {snp_filtered} "
            "--ts_filter_level 99.0 "
            "-mode SNP ").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                vcf=self.vcf,
                snp_tranches=self.snp_tranches,
                snp_filtered=self.snp_filtered,
                snp_recal=self.snp_recal)
       os.system("echo {0}".format(cmd))
       os.system("touch {0}".format(self.snp_filtered))

    def requires(self):
        return self.clone(VariantRecalibratorINDEL)

    def output(self):
        return luigi.LocalTarget(self.snp_filtered)

class ApplyRecalibrationINDEL(LocalSGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')

    def __init__(self, *args, **kwargs):
        super(ApplyRecalibrationINDEL, self).__init__(*args, **kwargs)
        self.snp_filtered = "{scratch}/{sample_name}/{sample_name}.snp.filtered.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_recal = "{scratch}/{sample_name}/{sample_name}.indel.recal".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_tranches = "{scratch}/{sample_name}/{sample_name}.indel.tranches".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_filtered = "{scratch}/{sample_name}/{sample_name}.snp.indel.filtered.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)

    def work(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T ApplyRecalibration "
            "-I {snp_filtered} "
            "-tranchesFile {indel_tranches} "
            "-recalFile {indel_recal} "
            "-o {indel_filtered} "
            "--ts_filter_level 99.0 "
            "-mode INDEL ").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                snp_filtered=self.snp_filtered,
                indel_recal=self.indel_recal,
                indel_tranches=self.indel_tranches,
                indel_filtered=self.indel_filtered)
       os.system("echo {0}".format(cmd))
       os.system("touch {0}".format(self.indel_filtered))

    def requires(self):
        return self.clone(ApplyRecalibrationSNP)

    def output(self):
        return luigi.LocalTarget(self.indel_filtered)

class VariantFiltrationSNP(LocalSGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')

    def __init__(self, *args, **kwargs):
        super(VariantFiltrationSNP, self).__init__(*args, **kwargs)
        self.vcf = "{scratch}/{sample_name}/{sample_name}.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.snp_filtered = "{scratch}/{sample_name}/{sample_name}.snp.filtered.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)

    def work(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T VariantFiltration "
            "-V {vcf} "
            "--filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" "
            "--filterName \"SNP_filter\" "
            "-o {snp_filtered} ").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                vcf=self.vcf,
                snp_filtered=self.snp_filtered)
       os.system("echo {0}".format(cmd))
       os.system("touch {0}".format(self.snp_filtered))

    def requires(self):
        return self.clone(HaplotypeCaller)

    def output(self):
        return luigi.LocalTarget(self.snp_filtered)

class VariantFiltrationINDEL(LocalSGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    java = luigi.Parameter(default=java,
        description = 'java version used')
    gatk = luigi.Parameter(default=gatk,
        description = 'gatk version used')
    max_mem = luigi.Parameter(default=max_mem,
        description = 'heap size for java in Gb')
    ref = luigi.Parameter(default=ref,
        description = 'reference genome location')

    def __init__(self, *args, **kwargs):
        super(VariantFiltrationINDEL, self).__init__(*args, **kwargs)
        self.snp_filtered = "{scratch}/{sample_name}/{sample_name}.snp.filtered.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_filtered = "{scratch}/{sample_name}/{sample_name}.snp.indel.filtered.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)

    def work(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T VariantFiltration "
            "-V {snp_filtered} "
            "--filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" "
            "--filterName \"INDEL_filter\" "
            "-o {indel_filtered}").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                snp_filtered=self.snp_filtered,
                indel_filtered=self.indel_filtered)
       os.system("echo {0}".format(cmd))
       os.system("touch {0}".format(self.indel_filtered))

    def requires(self):
        return self.clone(VariantFiltrationSNP)

    def output(self):
        return luigi.LocalTarget(self.indel_filtered)
