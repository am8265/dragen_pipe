import luigi
import os
import sys
import subprocess
from dragen_sample import dragen_sample
from luigi import configuration
from luigi.contrib.sge import SGEJobTask, LocalSGEJobTask

"""
Run samples through a luigized GATK pipeline after they have finished the
Dragen based alignment
"""

class config(luigi.Config):
    gatk = luigi.Parameter()
    java = luigi.Parameter()
    ref = luigi.Parameter()
    Mills1000g = luigi.Parameter()
    dbSNPLoc = luigi.Parameter()
    scratch = luigi.Parameter()
    max_mem = luigi.Parameter()
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    bam = luigi.Parameter()
    interval_list = luigi.Parameter()
    realn_bam = luigi.Parameter()
    recal_table = luigi.Parameter()
    vcf = luigi.Parameter()
    snp_vcf = luigi.Parameter()
    indel_vcf = luigi.Parameter()
    recal_snp = luigi.Parameter()
    recal_indel = luigi.Parameter()
    filter_snp = luigi.Parameter()
    filter_indel = luigi.Parameter()

    #scratch = getScratch()

    def getScratch():
        """get best scratch dir based on free space and # of pipeline jobs running"""
        pass

class RealignerTargetCreator(luigi.Task):
    """class for creating targets for indel realignment BAMs from Dragen"""

    def output(self):
        return luigi.LocalTarget(config().interval_list)

    def run(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T RealignerTargetCreator "
            "-I {bam} "
            "-o {interval_list} "
            "-known {Mills1000g} "
            "-known {dbSNPLoc} "
            "-nt 4 ").format(**config().__dict__)
        os.system("echo {0}".format(cmd))
        print luigi.LocalTarget(config().interval_list),config().interval_list

    def require():
        return luigi.LocalTarget(config().bam)

class IndelRealigner(luigi.Task):
    """class to create BAM with realigned BAMs"""

    def run(self):
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
            "-known {dbSNPLoc}").format(**config().__dict__)
        os.system("echo {0}".format(cmd))

    def require(self):
        return RealignerTargetCreator()

    def output(self):
        return luigi.LocalTarget(config().realn_bam)

class BaseRecalibrator(luigi.Target):

    def run(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T BaseRecalibrator "
            "-I {realn_bam) "
            "-nct 4 "
            "-o {recal_table} "
            "-known {Mills1000g} "
            "-known {dbSNPLoc}").format(**config().__dict__)

        os.system("echo {0}".format(cmd))

        rm_cmd = ['rm',config().bam]
        subprocess(rm_cmd)
    def require(self):
        return IndelRealigner

    def output(self):
        return luigi.LocalTarget(config().recal_table)

class PrintReads(luigi.Target):

    def run(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T PrintReads "
            "-I {realn_bam} "
            "-BQSR {recal_table} "
            "-o {realn_recal_bam} "
            "-nct 4").format(**config().__dict__)

    def require(self):
        return BaseRecalibrator

    def output(self):
        return luigi.LocalTarget(config().realn_recal_bam)

class HaplotypeCaller(luigi.Target):

    def run(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T HaplotypeCaller "
            "-I {realn_recal_bam} "
            "-o {vcf} "
            "-stand_call_conf 20 "
            "-stand_emit_conf 20 "
            "--dbsnp {dbSNP} "
            "-nct 4").format(**config().__dict__)

        rm_cmd = ['rm',config().realn_bam]
        subprocess(rm_cmd)

    def require(self):
        return BaseRecalibrator

    def output(self):
        return luigi.LocalTarget(config().realn_recal_bam)

class SelectVariantsSNP(luigi.Target):

    def run(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T SelectVariants "
            "-V {vcf} "
            "--selectType SNP "
            "-o {snp_vcf}").format(**config().__dict__)

    def require(self):
        return HaplotypeCaller

    def output(self):
        return luigi.LocalTarget(config().snp_vcf)

class SelectVariantsINDEL(luigi.Target):

    def run(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T SelectVariants "
            "-V {vcf} "
            "--selectType INDEL "
            "-o {indel_vcf}").format(**config().__dict__)

    def require(self):
        return HaplotypeCaller

    def output(self):
        return luigi.LocalTarget(config().indel_vcf)

class VariantRecalibratorSNP(luigi.Target):

    def run(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T VariantRecalibrator "
            "-I {snp_vcf} "
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
            "-tranchesFile {snp_tranches "
            "-rscriptFile {snp_rscript} "
            "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} "
            "-resource:omni,known=false,training=true,truth=false,prior=12.0 {omni} "
            "-resource:1000G,known=false,training=true,truth=false,prior=10.0 {1000g} "
            "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbSNP} "
            ).format(**config().__dict__)

    def require(self):
        return SelectVariantsSNP

    def output(self):
        return luigi.LocalTarget(config().snp_recal,config().snp_tranches,config().snp_rscript)

        ####double check resources version ####
        ####double check dbsnp version ####
class VariantRecalibratorINDEL(luigi.Target):

    def run(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T VariantRecalibrator "
            "-I {indel_vcf} "
            "-tranche 100.0 "
            "-tranche 99.9 "
            "-tranche 99.0 "
            "-tranche 90.0 "
            "-recalFile {indel_recal} "
            "-tranchesFile {indel_tranches "
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
            "-mode INDEL ").format(**config().__dict__)
    def require(self):
        return HaplotypeCaller

    def output(self):
        return luigi.LocalTarget(config().indel_recal,config().indel_tranches,config().indel_rscript)

class ApplyRecalibrationSNP(luigi.Target):

    def run(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T ApplyRecalibration "
            "-I {snp_vcf} "
            "-tranchesFile {snp_tranches} "
            "-recalFile {snp_recal} "
            "-o {snp_filtered} "
            "--ts_filter_level 99.0 "
            "-mode SNP ".format(**config().__dict__)

    def require(self):
        return VariantRecalibratorSNP

    def output(self):
        return luigi.LocalTarget(config().snp_vcf)

class ApplyRecalibrationINDEL(luigi.Target):

    def run(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T ApplyRecalibration "
            "-I {indel_vcf} "
            "-tranchesFile {indel_tranches} "
            "-recalFile {indel_recal} "
            "-o {indel_filtered} "
            "--ts_filter_level 99.0 "
            "-mode INDEL ".format(**config().__dict__)

    def require(self):
        return VariantRecalibratorINDEL

    def output(self):
        return luigi.LocalTarget(config().indel_vcf)

class VariantFiltrationSNP(luigi.Target):

    def run(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T VariantFiltration "
            "-V {snp_vcf} "
            "--filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" "
            "--filterName \"SNP_filter\" "
            "-o {snp_filtered} ").format(**config().__dict__)

    def require(self):
        return SelectVariantsSNP

    def output(self):
        return luigi.LocalTarget(config().snp_filtered)

class VariantFiltrationINDEL(luigi.Target):

    def run(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T VariantFiltration "
            "-V {indel_vcf} "
            "--filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" "
            "--filterName \"INDEL_filter\" "
            "-o {indel_filtered} "

    def require(self):
        return SelectVariantsINDEL

    def output(self):
        return luigi.LocalTarget(config().indel_filtered)

class IndelMatching(luigi.Target):

    def run(self):
        cmd = ""

    def require():
        if sample_type == 'genome':
            return ApplyRecalibrationINDEL
        else:
            return VariantFiltrationINDEL

    def output(self):
        return luigi.LocalTarget(config().indel_filter_matched)

"""
class CombineVariants(luigi.Target):

    def run(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
"""
