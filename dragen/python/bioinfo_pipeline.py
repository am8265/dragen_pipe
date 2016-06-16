#!/nfs/goldstein/software/python2.7.7/bin/python2.7
 
"""
Create a configuration file used by the Dragen system to processed IGM samples
based on their priority level in the Dragen pipeline queue
"""

import os
import shlex
import sys
import subprocess
import luigi
import MySQLdb
from luigi.contrib.sge import SGEJobTask, LocalSGEJobTask
from luigi.util import inherits
from dragen_sample import dragen_sample
import logging
import time


#Pipeline programs
gatk="/nfs/goldstein/software/GATK-3.6.0"
java="/nfs/goldstein/software/jdk1.8.0_05/bin/java"
bgzip="/nfs/goldstein/software/bin/bgzip"
tabix="/nfs/goldstein/software/bin/tabix"
picard="/nfs/goldstein/software/picard-tools-1.131/picard.jar"
bedtools="/nfs/goldstein/software/bedtools-2.25.0/bin/bedtools"

#GATK parameters
max_mem="15"
ref="/nfs/goldsteindata/refDB/HS_Build37/BWA_INDEX_hs37d5_BWAmem/hs37d5.fa"
hapmap="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/hapmap_3.3.b37.vcf"
omni="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/1000G_omni2.5.b37.vcf"
g1000="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/1000G_phase1.indels.b37.vcf"
Mills1000g="/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/Mills_and_1000G_gold_standard.indels.b37.vcf"
dbSNP="/nfs/goldsteindata/refDB/dbSNP/dbsnp_147.b37.vcf"
interval="/nfs/goldstein/software/dragen_pipe/dragen/conf/hs37d5.intervals"



class config(luigi.Config):
    """
    config class for instantiating parameters for this pipeline
    the values are read from luigi.cfg in the current folder 
    """

    '''
    ## GATK pipeline variables
    gatk = luigi.Parameter()
    #java = luigi.Parameter()
    max_mem = luigi.IntParameter()
    ref = luigi.Parameter()  
    scratch = luigi.Parameter()
    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    interval_list = luigi.Parameter()
    recal_table = luigi.Parameter()
    
    # BAMs
    bam = luigi.Parameter()
    realn_bam = luigi.Parameter()
    #realn_recal_bam = luigi.Parameter()

    
    # VQSR
    vcf = luigi.Parameter()
    snp_vcf = luigi.Parameter()
    indel_vcf = luigi.Parameter()
    
    # resources
    hapman=luigi.Parameter()
    omni=luigi.Parameter()
    g1000=luigi.Parameter()
    Mills1000g=luigi.Parameter()
    dbSNP=luigi.Parameter()
    '''

    ## Post GATK Variables,commented out redundant parameters
    java = luigi.Parameter()
    picard = luigi.Parameter()
    ref = luigi.Parameter()
    seqdict_file = luigi.Parameter()
    bed_file = luigi.Parameter()
    target_file = luigi.Parameter()
    create_targetfile = luigi.BooleanParameter()
    bait_file = luigi.Parameter()
    wgsinterval = luigi.BooleanParameter()
    scratch = luigi.Parameter(default='./')
    output_dir = luigi.Parameter()
    sample_name = luigi.Parameter()
    bam = luigi.Parameter()
    seq_type = luigi.Parameter()
    python_path = luigi.Parameter()
    relatedness_refs = luigi.Parameter()
    sampleped_loc = luigi.Parameter()
    relatedness_markers = luigi.Parameter()
    bedtools_loc = luigi.Parameter()
    pypy_loc = luigi.Parameter()
    coverage_binner_loc = luigi.Parameter()
    dbsnp = luigi.Parameter()


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
        self.bam = "{base_directory}/{sample_name}/{sample_name}.bam".format(
            base_directory=self.base_directory, sample_name=self.sample_name)
        self.interval_list = "{scratch}/{sample_name}/{sample_name}.interval_list".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
        os.system("mkdir -p {0}/scripts".format(self.scratch + self.sample_name))
        os.system("mkdir -p {0}/logs".format(self.scratch + self.sample_name))

        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T RealignerTargetCreator "
            "-L {interval} "
            "-I {bam} "
            "-o {interval_list} "
            "-known {Mills1000g} "
            "-known {dbSNP} "
            "-nt 4 ").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
                bam=self.bam,
                interval_list=self.interval_list,
                Mills1000g=Mills1000g,
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
        self.bam = "{base_directory}/{sample_name}/{sample_name}.bam".format(
            base_directory=self.base_directory, sample_name=self.sample_name)
        self.interval_list = "{scratch}/{sample_name}/{sample_name}.interval_list".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.realn_bam = "{scratch}/{sample_name}/{sample_name}.realn.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T IndelRealigner "
            "-L {interval} "
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
                interval=interval,
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
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk}/GenomeAnalysisTK.jar "
            "-R {ref} "
            "-T PrintReads "
            "-L {interval} "
            "-I {realn_bam} "
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
            "-L {interval} "
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
                interval=interval,
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
            "-o {snp_vcf}").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
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
            "-L {interval} "
            "-input {snp_filtered} "
            "-tranchesFile {indel_tranches} "
            "-recalFile {indel_recal} "
            "-o {indel_filtered} "
            "--ts_filter_level 99.0 "
            "-mode INDEL ").format(java=java,
                gatk=gatk,
                max_mem=max_mem,
                ref=ref,
                interval=interval,
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
        tabix_cmd = ("{tabix} {final_vcf_gz}").format(tabix=tabix,final_vcf_gz=self.final_vcf_gz)

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

class ArchiveSample:
    """Archive samples on Amplidata"""

    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter(default=interval)
    dragen_id = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    """
    db = MySQLdb.connect(db="sequenceDB", read_default_group="clientsequencedb",
            read_default_file="/nfs/goldstein/software/dragen/.my.cnf")
    curs = db.cursor()
    """

    def __init__(self, *args, **kwargs):
        super(ApplyRecalibrationSNP, self).__init__(*args, **kwargs)
        self.script = "{scratch}/{sample_name}/scripts/{class_name}.sh".format(
            scratch=self.scratch, sample_name=self.sample_name,class_name=self.__class__.__name__)
        self.recal_bam = "{scratch}/{sample_name}/{sample_name}.realn.recal.bam".format(
            scratch=self.scratch, sample_name=self.sample_name)
        self.bam = "{base_directory}/{sample_name}/{sample_name}.bam".format(
            base_directory=self.base_directory, sample_name=self.sample_name)
        self.final_vcf_gz = "{scratch}/{sample_name}/{sample_name}.analysisReady.vcf.gz".format(
            scratch=self.scratch, sample_name=self.sample_name)

    def work(self):
        cmd = ("rsync -a --partial --timeout=20000 -r "
              "{script_dir} {recal_bam} {final_vcf_gz}"
              ).format(recal_bam=recal_bam,
                      script_dir=script_dir,
                      final_vcf_gz=final_vcf_gz)

        with open(self.script,'w') as o:
            o.write(cmd + "\n")
            subprocess.check_call(shlex.split(cmd))

        # Original dragen BAM could be technically deleted earlier after the
        # realigned BAM has been created on scratch space but it is safer to 
        # delete after the final realigned, recalculated BAM has been archived
        # since our scratch space has failed in the past.
        rm_cmd = ['rm',self.bam,self.scratch_dir]
        rm_folder_cmd = ['rm','-rf',self.scratch_dir]
        #print rm_cmd
        #subprocess.call(rm_cmd)

    def requires(self):
        return self.clone(CombineVariants)

    def output(self):
        return luigi.LocalTarget(self.copy_complete)

################################################## END OF GATK PIPELINE ####################################################################################



##################################################################################################
##                                POST GATK TASKS                                               ##
##                                Added By : Raghav                                             ##
##################################################################################################
##   Notes :                                                                                    ##
##   1. Created a new config class, combining Josh and my variables                             ##
##   2. Working on integrating the two pipelines                                                ##   
##                                                                                              ##
##                                                                                              ##
##                                                                                              ##
##################################################################################################
##        Run as :                                                                              ##  
##        luigi --module luigized_coverage_statistics RunCvgMetrics  --local-scheduler          ##
##        Parameters are stored in luigi.cfg                                                    ##
##                                                                                              ##
##################################################################################################


class RootTask(luigi.WrapperTask):
    """
    Wrapper Task
    """

    ## Read in the bam file from the cmd line just to make sure
    ## to ensure uniqueness of the RootTask 
    bam = luigi.Parameter()
    
    
    def requires(self):
        self.output_file = "{0}/{1}.cvg.metrics.txt".format(self.scratch,self.sample_name)
        #self.output_file = config().scratch+config().sample_name+'.cvg.metrics.txt'
        return [ParsePicardOutputs(self.output_file,self.bam)]



class MyExtTask(luigi.ExternalTask):
    """
    Checks whether the file specified exists on disk
    """

    file_loc = luigi.Parameter()
    def output(self):
        return luigi.LocalTarget(self.file_loc)



class ParsePicardOutputs(SGEJobTask):
    """
    Simple Parser for Picard output
    """


    sample_name = luigi.Parameter()
    bam = luigi.Parameter()
    seq_type = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"

    def __init__(self,*args,**kwargs):
        """
        Initialize Class Variables
        """

        super(ParsePicardOutputs,self).__init__(*args,**kwargs)
        
        ## The input and output file to this task
        self.cvg_file = "{0}/{1}.cvg.metrics.txt".format(self.scratch,
                                                         self.sample)

        self.output_file = "{0}{1}.cvg.metrics.txt.transposed".format
        (config().scratch,self.sample)

        ## Define the command to be run 
        self.cmd = """cat %s | grep -v '^#'| head -3 | awk -f /home/rp2801/git/dragen_pipe/transpose.awk > %s """%(self.output_file,self.output_file)

    def requires(self):
        """ 
        Dependency is just the presence of the awk script
        and the completion of the RunCvgMetrics Task
        """
        
        return RunCvgMetrics(self.bam,self.sample_name,self.seq_type)

    def output(self):
        """
        Output from this task
        """

        return luigi.LocalTarget(self.output_file+'.transposed')

    def work(self):
        """
        Just a simple shell command is enough!
        """
        
        #os.system(cmd)
        os.system("""touch %s.transposed"""%(self.output_file))



class RunCvgMetrics(SGEJobTask):
    """ 
    A luigi task
    """
    
    
    bam = luigi.Parameter()
    sample_name = luigi.Parameter()
    seq_type = luigi.Parameter()
    
    
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"


    def __init__(self,*args,**kwargs):
        """
        Initialize Class Variables
        """

        super(RunCvgMetrics,self).__init__(*args,**kwargs)

        ## Initialize the output file
        self.output_file = "{0}/{1}.cvg.metrics.txt".format(
            config().scratch,self.sample_name)

    
        ## Initialize the targetfile name if it needs to be created
        if config().create_targetfile == True:
            self.bed_stem = utils.get_filename_from_fullpath(config().bed_file)
            self.target_file = "{0}/.{1}.list".format(config().scratch,self.bed_stem)
            #self.target_file = ((config().scratch) +
                           #(utils.get_filename_from_fullpath(config().bed_file)) +
                           #(".list"))
        else:
            self.target_file = config().target_file 

        ## Setup the command to be run 

        if self.seq_type.upper() == 'GENOME': ## If it is a genome sample
            if self.wgsinterval == True: ## Restrict wgs metrics to an interval
                cvg_cmd = ("%s -jar %s CollectWgsMetrics R=%s O=%s I=%s INTERVALS=%s MQ=20 Q=10"%(config().java,config().picard,config().reference_file,self.output_file,self.bam,self.target_file))
                
            else: ## Run wgs metrics across the genome
                cvg_cmd = ("%s -jar %s CollectWgsMetrics R=%s O=%s I=%s MQ=20 Q=10"%(config().java,config().picard,config().reference_file,self.output_file,self.bam))
              
                             
        else: ## Anything other than genome i.e. Exome or Custom Capture( maybe have an elif here to check further ?)
            cvg_cmd = ("%s -jar %s CalculateHsMetrics BI=%s TI=%s METRIC_ACCUMULATION_LEVEL=ALL_READS I=%s O=%s MQ=20 Q=10"%(config().java,config().picard,config().bait_file,self.target_file,self.bam,self.output_file))

        
    def output(self):
        """
        The output produced by this task 
        """
        return luigi.LocalTarget(self.output_file)
    

    def requires(self):
        """
        The dependency for this task is the CreateTargetFile task
        if the appropriate flag is specified by the user
        """

        if config().create_targetfile == True:
            yield [CreateTargetFile(),MyExtTask(self.bam)]
        else:
            yield MyExtTask(self.bam)
        
        
    def work(self):
        """
        Run Picard CalculateHsMetrics or WgsMetrics
        """
        ## Run the command
        os.system(self.cvg_cmd)
        #os.system("touch %s"%(self.output_file))


class CreateTargetFile(SGEJobTask):
    """ 
    A Luigi Task
    """
    

    def __init__(self,*args,**kwargs):
        super(CreateTargetFile,self).__init__(*args,**kwargs)

        self.bed_stem = utils.get_filename_from_fullpath(config().bed_file) + ".list"
        self.output_file = "{0}/{1}.genomecvg.bed".format(
                            config().scratch,self.bed_stem)

        self.bedtointerval_cmd = ("%s -jar %s BedToIntervalList I=%s SD=%s OUTPUT=%s"%(config().java,config().picard,config().bed_file,config().seqdict_file,self.output_file))    

    
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"



    def requires(self):
        """ 
        Dependcy for this task is the existence of the bedfile
        and the human build37 sequence dictionary file
        """
        
        return [MyExtTask(config().bed_file),MyExtTask(config().seqdict_file)]
    
    def output(self):
        """
        Returns the output from this task.
        In this case, a successful executation of this task will create a
        file on the local file system
        """      
        
        return luigi.LocalTarget(self.output_file)

    def work(self):
        """
        Run Picard BedToIntervalList 
        """
             
        #os.system("touch %s"%self.output_file)
        os.system(self.bedtointerval_cmd)

        
class CreatePed(SGEJobTask):
    """ 
    A luigi Task for creating ped files
    """
    
    ped_type = luigi.Parameter() ## This has one of two values :
                                 ## ethnicity or relatedness
    bam = luigi.Parameter()
    vcf = luigi.Parameter()
    sample_name = luigi.Parameter()
    family_id = luigi.Parameter(default="")
    seq_type = luigi.Parameter()
    
    
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"


    def __init__(self,*args,**kwargs):
        super(CreatePed,self).__init__(*args,**kwargs)

        self.output_file = config().scratch+self.sample_name+'.'+self.ped_type+'.ped'
        

        if self.ped_type == 'relatedness':
            self.sampletoped_cmd = ("%s %s --scratch %s --vcffile %s --bamfile %s" 
                               " --markers %s --refs %s --sample %s --fid %s"
                               " --seq_type %s --ped_type %s --output %s --pythondir /home/rp2801/git/dragen_pipe/"
                               %(config().python_path,config().sampleped_loc,
                                 config().scratch,self.vcf,
                                 self.bam,config().relatedness_markers,
                                 config().relatedness_refs,self.sample_name,
                                 self.family_id,self.seq_type,self.ped_type,
                                 config().scratch)
                               )
            
        elif self.ped_type == 'ethnicity':
            self.sampletoped_cmd = ("%s %s --scratch %s --vcffile %s --bamfile %s"
                               " --markers %s --refs %s --sample %s"
                               " --seq_type %s --ped_type %s --output %s" 
                               %(config().python_path,config().sampleped_loc,
                               config().scratch,self.vcf_file,
                               self.bam,config().relatedness_markers,
                               config().relatedness_refs,self.sample_name,
                               self.seq_type,self.ped_type,config().scratch)
                               )


    def requires(self):
        """
        The dependencies for this task are the existence of the relevant files
        """
        
        return [MyExtTask(self.bam),MyExtTask(self.vcf),MyExtTask(config().relatedness_refs)]


    def output(self):
        """
        Output from this task is a ped file named with correct type
        """
           
        return luigi.LocalTarget(self.output_file)


    def work(self):
        """
        To run this task, the sampletoped_fixed.py script is called
        """

            
        os.system(self.sampletoped_cmd)
        #os.system("%s sampletoped_fixed.py --scratch %s")

        
class CreateGenomeBed(SGEJobTask):
    """
    Use bedtools to create the input bed file for the coverage binning script
    """

    bam = luigi.Parameter()
    sample_name = luigi.Parameter()
    
    def __init__(self,*args,**kwargs):
        super(CreateGenomeBed,self).__init__(*args,**kwargs)

        self.genomecov_bed = "{0}{1}.genomecvg.bed".format(
                              config().scratch,self.sample)
        self.genomecov_cmd = "%s genomecov -bga -ibam %s > %s"%(
                              config().bedtools_loc,self.bam,
                              self.genomecov_bed)    
   
      
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"

    
    def requires(self):
        """
        The dependency is the presence of the bam file
        """
        
        return MyExtTask(self.bam)

    def output(self):
        """
        Output is the genomecov bed file
        """

        return luigi.LocalTarget(self.genomecov_bed)

    def work(self):
        """
        Run the bedtools cmd
        """

        os.system(self.genomecov_cmd)



class CoverageBinning(SGEJobTask):
    """
    Task to run the coverage binning script
    """


    bam_file = luigi.Parameter()
    sample = luigi.Parameter()
    block_size = luigi.Parameter()

    
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"
      
    
    def __init__(self,*args,**kwargs):
        super(CoverageBinning,self).__init__(*args,**kwargs)

        self.genomecov_bed = "{0}{1}.genomecvg.bed".format(
        config().scratch,self.sample)

        self.human_chromosomes = []
        self.human_chromosomes.extend(range(1, 23))
        self.human_chromosomes = [str(x) for x in self.human_chromosomes]
        self.human_chromosomes.extend(['X', 'Y','MT'])

        self.binning_cmd = "%s %s %s %s %s %s"%(config().pypy_loc,
                                     config().coverage_binner_loc,
                                     self.block_size,self.genomecov_bed,
                                     self.sample,config().scratch)
        
   
    def requires(self):
        """
        Dependency is the completion of the CreateGenomeBed Task
        """
       
        return [CreateGenomeBed(self.bam,self.sample_name)]
        
    
    def output(self):
        """
        Output from this task are 23 files for each chromosome
        """

        for chrom in self.human_chromosomes:
            yield [luigi.LocalTarget(config().scratch+self.sample+"_read_coverage_" +
                                     self.block_size + "_chr%s.txt"%chrom)]

    def work(self):
        """ Run the binning script
        """

        os.system(self.binning_cmd)


class AlignmentMetrics(LocalSGEJobTask):
    """
    Parse Alignment Metrics from dragen logs
    """
    

    bam = luigi.Parameter()
    sample_name = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"

    def __init__(self):
        super(AlignmentMetrics,self).__init__(*args,**kwargs)
        self.output_file = '{scratch}{samp}.alignment.metrics.txt'.format(scratch=config().scratch,samp=self.sample)
        self.cmd = "%s %s CollectAlignmentSummaryMetrics TMP_DIR=%s VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=%s INPUT=%s OUTPUT=%s"%(config().java,config().picard,config().scratch,config().ref,self.bam,self.output_file)               
    
    
    def exists(self):
        """
        Check whether the output file is present
        """

        return luigi.LocalTarget(self.output_file)

    def requires(self):
        """
        The dependencies for this task is simply the existence of the bam file
        from dragen with duplicates removed
        """
        
        return MyExtTask(self.bam) 

    def work(self):
        """
        Execute the command for this task 
        """

        os.system(self.cmd)


class DuplicateMetrics(SGEJobTask):
    """
    Parse Duplicate Metrics from dragen logs
    """
    
    dragen_log = luigi.Parameter()
    sample_name = luigi.Parameter()
    

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"

    
    def __init__(self):
        super(AlignmentMetrics,self).__init__(*args,**kwargs)
        self.output_file = "{0}/{1}.dups.tmp".format(config().scratch,self.sample_name)
        self.cmd = "grep 'duplicates marked' %s"%self.dragen_log
        

    def requires(self):
        """ 
        Dependencies for this task is the existence of the log file 
        """

        return MyExtTask(self.bam)


    def output(self):
        """
        """

        return luigi.LocalTarget(self.output_file)


    def work(self):
        """
        Execute this task
        """
        ## The regular expression to catch the percentage duplicates in the grep string
        catch_dup = re.compile('.*\((.*)%\).*')

        dragen_output = subprocess.check_output(self.cmd,shell=True)
        match = re.match(catch_dup,dragen_output)
        perc_duplicates = match.group(1)
        with open(self.output_file,'w') as OUT_HANDLE:
            print >> OUT_HANDLE,perc_duplicates



def EthnicityPipeline(LocalSGEJobTask):
    """
    Run the Ethnicity Pipeline
    """


class CheckAppendStatus(luigi.Target):
    """
    This is a target class for checking whether the database
    entry was correctly flagged
    """


    sample = luigi.Parameter()
    seq_type = luigi.Parameter()


    def __init__(self):
        super(CheckAppendStatus,self).__init__(*args,**kwargs)
        self.cmd = '''seqdb -e "use sequenceDB; SELECT to_append FROM table_name WHERE sample_name='%s',seq_type='%s' '''%(self.sample,self.seq_type)

    def exists(self):
        """
        Checks whether the database was updated
        """

        db_entry = subprocess.check_output(self.cmd)
        
        if db_entry != 'False':
            return False
        
        return True 


           
#@inherits(CreatePed) ## Try this method out too 
class CreateRelatednessPed(SGEJobTask):
    """
    Wrapper Task for creating a relatedness ped file and
    flagging the database entry to False
    """

    ped_type = luigi.Parameter(significant=True) ## This has one of two values :
                                 ## ethnicity or relatedness
    bam = luigi.Parameter(significant=True)
    #vcf_file = luigi.Parameter(significant=True)
    sample_name = luigi.Parameter(significant=True)
    family_id = luigi.Parameter(default="",significant=True)
    seq_type = luigi.Parameter(significant=True)
    master_ped = luigi.Parameter()


    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"
      
    
    def __init__(self,*args,**kwargs):
        super(CreateRelatednessPed,self).__init__(*args,**kwargs)       
        self.output_file = config().scratch+self.sample+'.'+self.ped_type+'.ped'        
        self.vcf = "{0}/{1}.analysisReady.vcf.gz".format(config().scratch,self.sample_name)
        self.cmd = '''seqdb -e "use sequenceDB; UPDATE table_name SET to_append=TRUE WHERE sample_name=%s, seq_type=%s" '''%(self.sample_name,self.seq_type)


    def requires(self):
        """
        The dependency for this task is the completion of the
        CreatePed Task with the appropriate signatures
        """
        
        return CreatePed(self.ped_type,self.bam,self.vcf,
                         self.sample_name,self.family_id,self.seq_type)

    
    def output(self):
        """
        The Output from this task 
        The database entry is flagged as False
        """
        
        return CheckAppendStatus(self.sample_name,self.seq_type)
        

    def work(self):
        """
        Run this task, a simple append
        should do it
        """

        os.system(self.cmd)

class AppendMasterPed(SGEJobTask):
    """
    Append the samples to the masterped
    """

    masterped = luigi.Parameter()
    threshold_file = luigi.Parameter()
    cryptic_threshold = luigi.Parameter()
    dbhost = luigi.Parameter()
    testmode = luigi.Parameter(default="2")
    
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"

    
    def __init__(self,*args,**kwargs):
        super(AppendMasterPed,self).__init__(*args,**kwargs)

        
    def requires(self):
        """
        Make sure relevant files are present
        """

        return [MyExtTask(self.masterped),MyExtTask(self.threshold_file)]
    
class RunRelatednessCheck(SGEJobTask):
    """
    This Task would be run nightly as a batch 
    """

    masterped = luigi.Parameter()
    #last_run = luigi.DateParameter()
        

class VariantCallingMetrics(SGEJobTask):
    """
    Run the picard tool for qc evaluation of the variant calls 
    """

    sample_name = luigi.Parameter()
    vcf = luigi.Parameter()
    

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"
      
    
    def __init__(self,*args,**kwargs):
        super(VariantCallingMetrics,self).__init__(*args,**kwargs)
        #sample_name = luigi.Parameter()
        self.output_file = "{0}/{1}.variant.metrics.txt".format(config().scratch,self.sample_name)
        self.cmd = "{java} -jar {picard} CollectVariantCallingMetrics INPUT={vcf} OUTPUT={out} DBSNP={db}".format(java=config().java,picard=config().picard,vcf=self.vcf,out=self.output_file,db=config().dbsnp)
        
       
    
    def requires(self):
        """
        The requirement for this task is the presence of the analysis ready 
        vcf file from the GATK pipeline
        """
    
        return MyExtTask(self.vcf)

    def output(self):
        """
        The result from this task is the creation of the metrics file
        """
        
        return luigi.LocalTarget(self.output_file)

    def work(self):
        """
        Run this task
        """
    
        
        os.system(self.cmd)
        #os.system("touch %s"%self.output_file) 


class UpdateDatabase(SGEJobTask):
    """
    Populate database with output files
    """

    output_file = luigi.Parameter()

    def __init__(self,*args,**kwargs):
        super(UpdateDatabase,self).__init__(*args,**kwargs)


    def requires(self):
        """
        The requirement for this task
        the outputfile from different tasks should
        be present
        """
        
        return [MyExtTask(),MyExtTask(),MyExtTask()]

    def output(self):
        """
        The output from this task
        """

    def work(self):
        """
        Run this task
        """

        

class GenderCheck(SGEJobTask):
    """
    Check gender using X,Y chromosome coverage
    """

    def __init__(self,*args,**kwargs):
        super(GenderCheck,self).__init__(*args,**kwargs)
        

    def requires(self):
        """
        The requirement for this task 
        """

        

    def output(self):
        """
        The output from this task
        """

    def work(self):
        """
        Run this task
        """
