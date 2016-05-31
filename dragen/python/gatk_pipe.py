#!/nfs/goldstein/software/python2.7.7/bin/python2.7

import os
import shlex
import sys
import subprocess
import luigi
from luigi.contrib.sge import SGEJobTask, LocalSGEJobTask
from dragen_sample import dragen_sample

"""
Run samples through a luigized GATK pipeline after they have finished the
Dragen based alignment
"""

gatk="/nfs/goldstein/software/GATK-3.5.0"
java="/nfs/goldstein/goldsteinlab/software/java/jdk1.7.0_03/bin/java"
bgzip="/nfs/goldstein/software/bin/bgzip"
tabix="/nfs/goldstein/software/bin/tabix"
picard="/nfs/goldstein/software/picard-tools-1.131/picard.jar"
max_mem="15"
ref="/nfs/goldsteindata/refDB/HS_Build37/BWA_INDEX_hs37d5_BWAmem/hs37d5.fa"
hapmap="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/hapmap_3.3.b37.vcf"
omni="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/1000G_omni2.5.b37.vcf"
g1000="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/1000G_phase1.indels.b37.vcf"
Mills1000g="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/Mills_and_1000G_gold_standard.indels.b37.vcf"
dbSNP="/nfs/goldsteindata/refDB/dbSNP/dbsnp_147.b37.vcf"
chr_list="/nfs/goldstein/software/dragen_pipe/dragen/conf/chr_list.bed"

class RootTask(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    chr_list = luigi.Parameter()

    def requires(self):
        if self.sample_type == 'genome':
            return self.clone(IndelRealigner)
    def work(self):
        os.system('touch /home/jb3816/{0}.complete'.format(self.sample_name))
    def output(self):
        return luigi.LocalTarget('/home/jb3816/{0}.complete'.format(self.sample_name))

class RealignerTargetCreator(SGEJobTask):
    """class for creating targets for indel realignment BAMs from Dragen"""

    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

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
        self.bam = "{base_directory}/{sample_name}/{sample_name}.bam".format(
            base_directory=self.base_directory, sample_name=self.sample_name)
        self.interval_list = "{scratch}/{sample_name}/{sample_name}.interval_list".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/RealignerTargetCreator.sh".format(
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
                chr_list=chr_list,
                dbSNP=dbSNP)
        if not os.path.isdir(os.path.dirname(self.script)):
            os.makedirs(os.path.dirname(self.script))
        with open(self.script,'w') as o:
            o.write(cmd + "\n")
        subprocess.check_call(shlex.split(cmd))

    def output(self):
        yield luigi.LocalTarget(self.interval_list)

class IndelRealigner(SGEJobTask):
    """class to create BAM with realigned BAMs"""
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

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
        self.bam = "{base_directory}/{sample_name}/{sample_name}.bam".format(
            base_directory=self.base_directory, sample_name=self.sample_name)
        self.interval_list = "{scratch}/{sample_name}/{sample_name}.interval_list".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.realn_bam = "{scratch}/{sample_name}/{sample_name}.realn.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/IndelRealigner.sh".format(
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
        self.script = "{scratch}/{sample_name}/scripts/script.sh".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)


    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T BaseRecalibrator "
            "-I {realn_bam} "
            "-nct 4 "
            "-o {recal_table} "
            "-L {capture_kit_bed} "
            "-knownSites {Mills1000g} "
            "-knownSites {dbSNP}").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                realn_bam=self.realn_bam,
                recal_table=self.recal_table,
                Mills1000g=Mills1000g,
                capture_kit_bed=self.capture_kit_bed,
                dbSNP=dbSNP)
        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))

        #os.system("echo {0}".format(cmd))
        #os.system("touch {0}".format(self.recal_table))

        rm_cmd = ['rm',self.bam]
        #print rm_cmd
        #subprocess.call(rm_cmd)

    def requires(self):
        return self.clone(IndelRealigner)

    def output(self):
        return luigi.LocalTarget(self.recal_table)

class PrintReads(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

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

        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))

    def requires(self):
        return self.clone(BaseRecalibrator)

    def output(self):
        return luigi.LocalTarget(self.recal_bam)

class HaplotypeCaller(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

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
    n_cpu = 4
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(HaplotypeCaller, self).__init__(*args, **kwargs)
        self.recal_table = "{scratch}/{sample_name}/{sample_name}.recal_table".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.recal_bam = "{scratch}/{sample_name}/{sample_name}.realn.recal.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.gvcf = "{scratch}/{sample_name}/{sample_name}.g.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T HaplotypeCaller "
            "-I {recal_bam} "
            "-o {gvcf} "
            "-stand_call_conf 20 "
            "-stand_emit_conf 20 "
            "--emitRefConfidence GVCF "
            "--variant_index_type LINEAR "
            "--variant_index_parameter 128000 "
            "--dbsnp {dbSNP} "
            "-nct 4").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                recal_bam=self.recal_bam,
                gvcf=self.gvcf,
                dbSNP=dbSNP)

        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))

       #rm_cmd = ['rm',self.recal_bam]
       #print rm_cmd
       #subprocess.call(rm_cmd)

    def requires(self):
        return self.clone(PrintReads)

    def output(self):
        return luigi.LocalTarget(self.gvcf)

class GenotypeGVCFs(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

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
        self.recal_bam = "{scratch}/{sample_name}/{sample_name}.realn.recal.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.gvcf = "{scratch}/{sample_name}/{sample_name}.g.vcf".format(
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
            "-o {vcf} "
            "-stand_call_conf 20 "
            "-stand_emit_conf 20 "
            "-V {gvcf}").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                gvcf=self.gvcf,
                vcf=self.vcf)

        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))

       #rm_cmd = ['rm',self.recal_bam]
       #print rm_cmd
       #subprocess.call(rm_cmd)

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
            "-V {vcf} "
            "-selectType SNP "
            "-o {snp_vcf}").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
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
            "-V {vcf} "
            "-selectType INDEL "
            "-o {snp_vcf}").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                snp_vcf=self.indel_vcf,
                vcf=self.vcf)

        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))


    def requires(self):
        return self.clone(SelectVariantsSNP)

    def output(self):
        return luigi.LocalTarget(self.indel_vcf)


class VariantRecalibratorSNP(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

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
        if self.sample_type == 'genome':
            return self.clone(SelectVariantsSNP)
        elif self.sample_type == 'exome':
            return self.clone(SelectVariantsSNP)
        else:
            raise Exception, "Sample type: %s not supported in this module" % self.sample_type


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
        self.indel_tranches = "{scratch}/{sample_name}/{sample_name}.indel.tranches".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T VariantRecalibrator "
            "-I {indel_vcf} "
            "-an QD "
            "-an DP "
            "-an FS "
            "-an SOR "
            "-an MQRankSum "
            "-an ReadPosRankSum "
            "-an InbreedingCoeff "
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
                indel_vcf=self.indel_vcf,
                indel_recal=self.indel_recal,
                indel_tranches=self.indel_tranches,
                indel_rscript=self.indel_rscript,
                Mills1000g=Mills1000g)
        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))

    def requires(self):
        if self.sample_type == 'genome':
            return self.clone(VariantRecalibratorSNP)
        else:
            raise Exception, "Sample type: %s not supported in this module" % self.sample_type

    def output(self):
        return luigi.LocalTarget(self.indel_recal)

class ApplyRecalibrationSNP(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

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
            "-input {vcf} "
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
        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))

    def requires(self):
        if self.sample_type == 'genome':
            return self.clone(VariantRecalibratorINDEL)
        elif self.sample_type == 'exome':
            return self.clone(VariantRecalibratorSNP)
        else:
            raise Exception, "Sample type: %s not supported in this module" % self.sample_type

    def output(self):
        return luigi.LocalTarget(self.snp_filtered)

class ApplyRecalibrationINDEL(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

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
        self.snp_filtered = "{scratch}/{sample_name}/{sample_name}.snp.filtered.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_recal = "{scratch}/{sample_name}/{sample_name}.indel.recal".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_tranches = "{scratch}/{sample_name}/{sample_name}.indel.tranches".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.indel_filtered = "{scratch}/{sample_name}/{sample_name}.snp.indel.filtered.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
       cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T ApplyRecalibration "
            "-input {snp_filtered} "
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
       with open(self.script,'w') as o:
           o.write(cmd + "\n")
           subprocess.check_call(shlex.split(cmd))

    def requires(self):
        if self.sample_type == 'genome':
            return self.clone(ApplyRecalibrationSNP)
        else:
            raise Exception, "Sample type: %s not supported in this module" % self.sample_type

    def output(self):
        return luigi.LocalTarget(self.indel_filtered)

class VariantFiltrationSNP(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

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
            "-V {snp_vcf} "
            "--filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" "
            "--filterName \"SNP_filter\" "
            "-o {snp_filtered} ").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                snp_vcf=self.snp_vcf,
                snp_filtered=self.snp_filtered)
       with open(self.script,'w') as o:
           o.write(cmd + "\n")
           subprocess.check_call(shlex.split(cmd))

    def requires(self):
        if self.sample_type == 'custom_capture':
            return self.clone(SelectVariantsSNP)
        else:
            raise Exception, "Sample type: %s not supported in this module" % self.sample_type

    def output(self):
        return luigi.LocalTarget(self.snp_filtered)

class VariantFiltrationINDEL(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

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
            "-V {indel_vcf} "
            "--filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" "
            "--filterName \"INDEL_filter\" "
            "-o {indel_filtered}").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                indel_vcf=self.indel_vcf,
                indel_filtered=self.indel_filtered)
       with open(self.script,'w') as o:
           o.write(cmd + "\n")
           subprocess.check_call(shlex.split(cmd))

    def requires(self):
        if self.sample_type == 'exome':
            return self.clone(ApplyRecalibrationSNP),self.clone(SelectVariantsINDEL)
        elif self.sample_type == 'custom_capture':
            return self.clone(VariantFiltrationINDEL)
        else:
            raise Exception, "Sample type: %s not supported in this module" % self.sample_type
    def output(self):
        return luigi.LocalTarget(self.indel_filtered)


class CombineVariants(SGEJobTask):
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

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
        self.final_vcf_gz = "{scratch}/{sample_name}/{sample_name}.analysisReady.vcf.gz".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.tmp_vcf = "{scratch}/{sample_name}/{sample_name}.tmp.vcf".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
        filter_flag = 0
        with open(self.tmp_vcf,'w') as vcf_out:
            with open(self.snp_filtered) as header:
                for line in header.readlines():
                    if line[0] == '#':
                        if line[0:8] == '##FILTER' and filter_flag == 0 :
                            filter_flag = 1
                            #SNP specific filters
                            if self.sample_type == 'exome' or self.sample_type =='genome':
                                vcf_out.write('##FILTER=<ID=VQSRTrancheSNP99.00to99.90,Description="Truth sensitivity tranche level for SNP model\n')
                                vcf_out.write('##FILTER=<ID=VQSRTrancheSNP99.90to100.00+,Description="Truth sensitivity tranche level for SNP model\n')
                                vcf_out.write('##FILTER=<ID=VQSRTrancheSNP99.90to100.00,Description="Truth sensitivity tranche level for SNP model\n')
                            else:
                                vcf_out.write('##FILTER=<ID=SNP_filter,Description="QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0">\n')

                            #Indel specific filters
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

        pre_cmd = ("/bin/bash -c cat <(head -1000 {vcf} | grep ^#) "
            "<(grep -v ^# {snp_filtered}) "
            "<(grep -v ^# {indel_filtered}) >> "
            "{tmp_vcf} "
            ).format(snp_filtered=self.snp_filtered,
                indel_filtered=self.indel_filtered,
                tmp_vcf=self.tmp_vcf,
                vcf=self.vcf)

        sort_cmd = ("{java} -jar {picard} "
            "SortVcf "
            "I={tmp_vcf} "
            "O={final_vcf}").format(java=java,
                picard=picard,
                tmp_vcf=self.tmp_vcf,
                final_vcf=self.final_vcf)

        bgzip_cmd = ("{bgzip} {final_vcf}").format(bgzip=bgzip,final_vcf=self.final_vcf)
        tabix_cmd = ("{tabix} {final_vcf}.gz").format(tabix=tabix,final_vcf=self.final_vcf)

        #rm_cmd = ('').format()
        with open(self.script,'w') as o:
            o.write(pre_cmd + "\n")
            o.write(sort_cmd + "\n")
            o.write(bgzip_cmd + "\n")
            o.write(tabix_cmd + "\n")
            #o.write(rm_cmd + "\n")

        subprocess.check_call(shlex.split(sort_cmd))
        subprocess.check_call(shlex.split(bgzip_cmd))
        subprocess.check_call(shlex.split(tabix_cmd))

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
        return luigi.LocalTarget(self.final_vcf_gz)
