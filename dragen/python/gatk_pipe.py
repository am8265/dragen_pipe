#!/nfs/goldstein/software/python2.7.7/bin/python2.7

import getpass
import os
import shlex
import sys
import subprocess
import luigi
import MySQLdb
from luigi.contrib.sge import SGEJobTask
from db_statements import *
from dragen_globals import *
from dragen_sample import dragen_sample

"""Run samples through a luigized GATK pipeline after finishing the
Dragen based alignment"""

def getDBIDMaxPrepID(pseudo_prepid):
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

def getUserID():
    userName = getpass.getuser()
    db = get_connection("seqdb")
    try:
        cur = db.cursor()
        cur.execute(GET_USERID.format(
            userName=userName))
        userid = cur.fetchone()
    finally:
        if db.open:
            db.close()
    return userid[0]

def getCaptureKitBed(pseudo_prepid):
    db = get_connection("seqdb")
    try:
        cur = db.cursor()
        cur.execute(GET_CAPTURE_KIT_BED.format(
            pseudo_prepid=pseudo_prepid))
        userid = cur.fetchone()
    finally:
        if db.open:
            db.close()
    return capture_kit_bed

class config(luigi.Config):
    """config class for instantiating parameters for this pipeline
    the values are read from luigi.cfg in the current folder"""
    #programs
    bedtools = luigi.Parameter()
    bgzip = luigi.Parameter()
    coverage_binner_loc = luigi.Parameter()
    gatk = luigi.Parameter()
    java = luigi.Parameter()
    picard = luigi.Parameter()
    pypy_loc = luigi.Parameter()
    snpEff = luigi.Parameter()
    tabix = luigi.Parameter()
    #files
    bait_file = luigi.Parameter()
    bait_file_X = luigi.Parameter()
    bait_file_Y = luigi.Parameter()
    ccds_bed_file = luigi.Parameter()
    relatedness_markers = luigi.Parameter()
    relatedness_refs = luigi.Parameter()
    sampleped_loc = luigi.Parameter()
    snpEff_cfg = luigi.Parameter()
    snpEff_interval = luigi.Parameter()
    target_file = luigi.Parameter()
    #gatk resources
    ref = luigi.Parameter()
    seqdict_file = luigi.Parameter()
    interval = luigi.Parameter()
    g1000 = luigi.Parameter()
    hapmap = luigi.Parameter()
    Mills1000g = luigi.Parameter()
    omni = luigi.Parameter()
    dbSNP = luigi.Parameter()
    #variables
    create_targetfile = luigi.BooleanParameter()
    max_mem = luigi.IntParameter()
    block_size = luigi.Parameter(description='The block size over which to do the binning')
    #locations
    python_path = luigi.Parameter()
    scratch = luigi.Parameter()
    base = luigi.Parameter()

class db(luigi.Config):
    """Database config variable will be read from the
    db section of the config file"""
    cnf = luigi.Parameter()
    seqdb_group = luigi.Parameter()
    dragen_group = luigi.Parameter()

class CopyBam(SGEJobTask):
    """class for copying dragen aligned bam to a scratch location"""
    sample_name = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"
    #job_name = self.sample_name + '.' + self.__class__.__name__ 

    def __init__(self, *args, **kwargs):
        super(CopyBam, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.base_dir = "{0}/{1}/{2}.{3}".format(
            config().base,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.bam = "{0}/{1}.{2}.bam".format(
            self.base_dir,self.sample_name,self.pseudo_prepid)
        self.bam_index = "{0}/{1}.{2}.bam.bai".format(
            self.base_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

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

            os.system("mkdir -p {0}/scripts".format(self.scratch_dir))
            os.system("mkdir -p {0}/logs".format(self.scratch_dir))
            cmd = ("rsync -a --timeout=20000 -r {bam} {bam_index} {scratch_dir}/"
                  ).format(bam=self.bam,
                           bam_index=self.bam_index,
                           scratch_dir=self.scratch_dir)
            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
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

class RealignerTargetCreator(SGEJobTask):
    """class for creating targets for indel realignment BAMs from Dragen"""
    sample_name = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()

    n_cpu = 4
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(RealignerTargetCreator, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.interval_list = "{0}/{1}.{2}.interval_list".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.scratch_bam = "{0}/{1}.{2}.bam".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

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
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T RealignerTargetCreator "
            "-L {interval} "
            "-I {scratch_bam} "
            "-o {interval_list} "
            "-known {Mills1000g} "
            "-known {dbSNP} "
            "--log_to_file {log_file} "
            "-nt 4 ").format(java=config().java,
                gatk=config().gatk,
                max_mem=config().max_mem,
                ref=config().ref,
                interval=config().interval,
                scratch_bam=self.scratch_bam,
                interval_list=self.interval_list,
                Mills1000g=config().Mills1000g,
                dbSNP=config().dbSNP,
                log_file=self.log_file)

        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(CopyBam)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)


class IndelRealigner(SGEJobTask):
    """class to create BAM with realigned BAMs"""
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(IndelRealigner, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.scratch_bam = "{0}/{1}.{2}.bam".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.interval_list = "{0}/{1}.{2}.interval_list".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.realn_bam = "{0}/{1}.{2}.realn.bam".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)

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
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T IndelRealigner "
            "-L {interval} "
            "-I {scratch_bam} "
            "-o {realn_bam} "
            "--log_to_file {log_file} "
            "-targetIntervals {interval_list} "
            "-maxReads 10000000 "
            "-maxInMemory 450000 "
            "-known {Mills1000g} "
            "-known {dbSNP}").format(java=config().java,
                gatk=config().gatk,
                max_mem=config().max_mem,
                ref=config().ref,
                interval=config().interval,
                scratch_bam=self.scratch_bam,
                realn_bam=self.realn_bam,
                interval_list=self.interval_list,
                Mills1000g=config().Mills1000g,
                log_file=self.log_file,
                dbSNP=config().dbSNP)
        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()

        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(RealignerTargetCreator)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class BaseRecalibrator(SGEJobTask):
    """class to create a recalibration table with realigned BAMs"""
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 4
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(BaseRecalibrator, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.realn_bam = "{0}/{1}.{2}.realn.bam".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.recal_table = "{0}/{1}.{2}.recal_table".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)

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
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T BaseRecalibrator "
            "-L {interval} "
            "-I {realn_bam} "
            "-o {recal_table} "
            "-L {capture_kit_bed} "
            "--log_to_file {log_file} "
            "-knownSites {Mills1000g} "
            "-knownSites {dbSNP} "
            "-nct 4").format(java=config().java,
                gatk=config().gatk,
                max_mem=config().max_mem,
                ref=config().ref,
                interval=config().interval,
                realn_bam=self.realn_bam,
                recal_table=self.recal_table,
                Mills1000g=config().Mills1000g,
                log_file=self.log_file,
                capture_kit_bed=self.capture_kit_bed,
                dbSNP=config().dbSNP)

        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()

        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(IndelRealigner)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class PrintReads(SGEJobTask):
    """class to create BAM with realigned and recalculated BAMs"""
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 4
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(PrintReads, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.realn_bam = "{0}/{1}.{2}.realn.bam".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.recal_table = "{0}/{1}.{2}.recal_table".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.recal_bam = "{0}/{1}.{2}.realn.recal.bam".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)

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
        # --disable_indel_quals are necessary to remove BI and BD tags in the bam file
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T PrintReads "
            "-L {interval} "
            "-I {realn_bam} "
            "--log_to_file {log_file} "
            "--disable_indel_quals "
            "-BQSR {recal_table} "
            "-o {recal_bam} "
            "-nct 4").format(java=config().java,
                gatk=config().gatk,
                max_mem=config().max_mem,
                ref=config().ref,
                interval=config().interval,
                log_file=self.log_file,
                realn_bam=self.realn_bam,
                recal_table=self.recal_table,
                recal_bam=self.recal_bam)

        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(BaseRecalibrator)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)


class HaplotypeCaller(SGEJobTask):
    """class to create gVCFs"""
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 4
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(HaplotypeCaller, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.realn_bam = "{0}/{1}.{2}.realn.bam".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.recal_table = "{0}/{1}.{2}.recal_table".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.scratch_bam = "{0}/{1}.{2}.bam".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.recal_bam = "{0}/{1}.{2}.realn.recal.bam".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.gvcf = "{0}/{1}.{2}.g.vcf.gz".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.gvcf_index = "{0}/{1}.{2}.g.vcf.gz.tbi".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

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
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T HaplotypeCaller "
            "-L {interval} "
            "-I {recal_bam} "
            "-o {gvcf} "
            "--log_to_file {log_file} "
            "-stand_call_conf 20 "
            "-stand_emit_conf 20 "
            "--emitRefConfidence GVCF "
            "-GQB 5 -GQB 15 -GQB 20 -GQB 60 "
            "--variant_index_type LINEAR "
            "--variant_index_parameter 128000 "
            "--dbsnp {dbSNP} "
            "-nct 4").format(java=config().java,
                gatk=config().gatk,
                max_mem=config().max_mem,
                ref=config().ref,
                interval=config().interval,
                log_file=self.log_file,
                recal_bam=self.recal_bam,
                gvcf=self.gvcf,
                dbSNP=config().dbSNP)

        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
        finally:
            if db.open:
                db.close()

        rm_cmd = ['rm',self.realn_bam,self.scratch_bam]
        subprocess.call(rm_cmd)

    def requires(self):
        return self.clone(PrintReads)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class GenotypeGVCFs(SGEJobTask):
    """class to perfrom variant calling from gVCFs"""
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(GenotypeGVCFs, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.recal_table = "{0}/{1}.{2}.recal_table".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.gvcf = "{0}/{1}.{2}.g.vcf.gz".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.vcf = "{0}/{1}.{2}.raw.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

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
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T GenotypeGVCFs "
            "-L {interval} "
            "-o {vcf} "
            "--log_to_file {log_file} "
            "-stand_call_conf 20 "
            "-stand_emit_conf 20 "
            "-V {gvcf}").format(java=config().java,
                gatk=config().gatk,
                max_mem=config().max_mem,
                log_file=self.log_file,
                ref=config().ref,
                interval=config().interval,
                gvcf=self.gvcf,
                vcf=self.vcf)

        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(HaplotypeCaller)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class SelectVariantsSNP(SGEJobTask):
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(SelectVariantsSNP, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.recal_table = "{0}/{1}.{2}.recal_table".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.vcf = "{0}/{1}.{2}.raw.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.snp_vcf = "{0}/{1}.{2}.snp.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

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
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T SelectVariants "
            "-L {interval} "
            "-V {vcf} "
            "--log_to_file {log_file} "
            "-selectType SNP "
            "-o {snp_vcf}").format(java=config().java,
                gatk=config().gatk,
                max_mem=config().max_mem,
                log_file=self.log_file,
                ref=config().ref,
                interval=config().interval,
                snp_vcf=self.snp_vcf,
                vcf=self.vcf)

        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()

        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(GenotypeGVCFs)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class SelectVariantsINDEL(SGEJobTask):
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(SelectVariantsINDEL, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.recal_table = "{0}/{1}.{2}.recal_table".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.vcf = "{0}/{1}.{2}.raw.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.indel_vcf = "{0}/{1}.{2}.indel.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

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
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T SelectVariants "
            "-L {interval} "
            "-V {vcf} "
            "--log_to_file {log_file} "
            "-selectType INDEL "
            "-o {indel_vcf}").format(java=config().java,
                gatk=config().gatk,
                max_mem=config().max_mem,
                log_file=self.log_file,
                ref=config().ref,
                interval=config().interval,
                indel_vcf=self.indel_vcf,
                vcf=self.vcf)

        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()

        finally:
            if db.open:
                db.close()

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
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)


class VariantRecalibratorSNP(SGEJobTask):
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(VariantRecalibratorSNP, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.snp_vcf = "{0}/{1}.{2}.snp.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.snp_recal = "{0}/{1}.{2}.snp.recal".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.snp_tranches = "{0}/{1}.{2}.snp.tranches".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.snp_rscript = "{0}/{1}.{2}.snp.rscript".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

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
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T VariantRecalibrator "
            "-L {interval} "
            "--log_to_file {log_file} "
            "--input {snp_vcf} "
            "-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum "
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
            ).format(java=config().java,
                gatk=config().gatk,
                max_mem=config().max_mem,
                log_file=self.log_file,
                ref=config().ref,
                interval=config().interval,
                snp_vcf=self.snp_vcf,
                snp_recal=self.snp_recal,
                snp_tranches=self.snp_tranches,
                snp_rscript=self.snp_rscript,
                dbSNP=config().dbSNP,
                hapmap=config().hapmap,
                omni=config().omni,
                g1000=config().g1000)

        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()

        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(SelectVariantsSNP)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

        ####double check dbsnp version ####

class VariantRecalibratorINDEL(SGEJobTask):
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(VariantRecalibratorINDEL, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.indel_vcf = "{0}/{1}.{2}.indel.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.indel_recal = "{0}/{1}.{2}.indel.recal".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.indel_rscript = "{0}/{1}.{2}.indel.rscript".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.indel_tranches = "{0}/{1}.{2}.indel.tranches".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

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
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T VariantRecalibrator "
            "-L {interval} "
            "--log_to_file {log_file} "
            "--input {indel_vcf} "
            "-an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum "
            "-mode INDEL "
            "--maxGaussians 4 "
            "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 "
            "-recalFile {indel_recal} "
            "-tranchesFile {indel_tranches} "
            "-resource:mills,known=true,training=true,truth=true,prior=12.0 {Mills1000g} "
            ).format(java=config().java,
                gatk=config().gatk,
                max_mem=config().max_mem,
                log_file=self.log_file,
                ref=config().ref,
                interval=config().interval,
                indel_vcf=self.indel_vcf,
                indel_recal=self.indel_recal,
                indel_tranches=self.indel_tranches,
                indel_rscript=self.indel_rscript,
                Mills1000g=config().Mills1000g)

        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()

        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(SelectVariantsINDEL)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class ApplyRecalibrationSNP(SGEJobTask):
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(ApplyRecalibrationSNP, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.vcf = "{0}/{1}.{2}.snp.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.snp_recal = "{0}/{1}.{2}.snp.recal".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.snp_tranches = "{0}/{1}.{2}.snp.tranches".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.snp_filtered = "{0}/{1}.{2}.snp.filtered.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

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
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T ApplyRecalibration "
            "-L {interval} "
            "--log_to_file {log_file} "
            "-input {vcf} "
            "-tranchesFile {snp_tranches} "
            "-recalFile {snp_recal} "
            "-o {snp_filtered} "
            "--ts_filter_level 99.0 "
            "-mode SNP ").format(java=config().java,
                gatk=config().gatk,
                max_mem=config().max_mem,
                ref=config().ref,
                log_file=self.log_file,
                interval=config().interval,
                vcf=self.vcf,
                snp_tranches=self.snp_tranches,
                snp_filtered=self.snp_filtered,
                snp_recal=self.snp_recal)

        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()

        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(VariantRecalibratorSNP)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class ApplyRecalibrationINDEL(SGEJobTask):

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(ApplyRecalibrationINDEL, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.indel_vcf = "{0}/{1}.{2}.indel.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.indel_recal = "{0}/{1}.{2}.indel.recal".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.indel_tranches = "{0}/{1}.{2}.indel.tranches".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.indel_filtered = "{0}/{1}.{2}.indel.filtered.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

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
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T ApplyRecalibration "
            "-L {interval} "
            "--log_to_file {log_file} "
            "-input {indel_vcf} "
            "-tranchesFile {indel_tranches} "
            "-recalFile {indel_recal} "
            "-o {indel_filtered} "
            "--ts_filter_level 99.0 "
            "-mode INDEL ").format(java=config().java,
                gatk=config().gatk,
                max_mem=config().max_mem,
                ref=config().ref,
                log_file=self.log_file,
                interval=config().interval,
                indel_vcf=self.indel_vcf,
                indel_recal=self.indel_recal,
                indel_tranches=self.indel_tranches,
                indel_filtered=self.indel_filtered)

        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()

        finally:
            if db.open:
                db.close()

    def requires(self):
      return self.clone(VariantRecalibratorINDEL)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class VariantFiltrationSNP(SGEJobTask):

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(VariantFiltrationSNP, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.snp_vcf = "{0}/{1}.{2}.snp.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.snp_filtered = "{0}/{1}.{2}.snp.filtered.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

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
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T VariantFiltration "
            "-L {interval} "
            "-V {snp_vcf} "
            "--filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" "
            "--filterName \"SNP_filter\" "
            "--log_to_file {log_file} "
            "-o {snp_filtered} ").format(java=config().java,
                gatk=config().gatk,
                max_mem=config().max_mem,
                ref=config().ref,
                log_file=self.log_file,
                interval=config().interval,
                snp_vcf=self.snp_vcf,
                snp_filtered=self.snp_filtered)

        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()

        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(SelectVariantsSNP)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class VariantFiltrationINDEL(SGEJobTask):

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(VariantFiltrationINDEL, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.indel_vcf = "{0}/{1}.{2}.indel.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.indel_filtered = "{0}/{1}.{2}.indel.filtered.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

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
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T VariantFiltration "
            "-L {interval} "
            "-V {indel_vcf} "
            "--filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" "
            "--filterName \"INDEL_filter\" "
            "--log_to_file {log_file} "
            "-o {indel_filtered}").format(java=config().java,
                gatk=config().gatk,
                max_mem=config().max_mem,
                ref=config().ref,
                log_file=self.log_file,
                interval=config().interval,
                indel_vcf=self.indel_vcf,
                indel_filtered=self.indel_filtered)

        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()

        finally:
            if db.open:
                db.close()

    def requires(self):
      return self.clone(SelectVariantsINDEL)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class CombineVariants(SGEJobTask):

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(CombineVariants, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.snp_filtered = "{0}/{1}.{2}.snp.filtered.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.vcf = "{0}/{1}.{2}.raw.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.indel_filtered = "{0}/{1}.{2}.indel.filtered.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.final_vcf = "{0}/{1}.{2}.analysisReady.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.tmp_vcf = "{0}/{1}.{2}.tmp.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

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
        """Merges SNP and INDEL vcfs.  Using the filtered SNP vcf header as a
        base, variant type/sample type specific ##FILTERs are added to the header.
        Additionally the AnnoDBID annotation is added.  After finishing reading
        through the header the SNP vcf is read again with only variants
        being outputted.  The same happens with the INDEL vcf.  The resulting
        vcf is then sorted producing the analysisReady.vcf.
        """
        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

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
                "O={final_vcf}").format(java=config().java,
                    picard=config().picard,
                    tmp_vcf=self.tmp_vcf,
                    final_vcf=self.final_vcf)
            with open(self.script,'w') as o:
                o.write(sort_cmd + "\n")
            subprocess.check_call(shlex.split(sort_cmd))

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()

        finally:
            if db.open:
                db.close()

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
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class AnnotateVCF(SGEJobTask):

    sample_name = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(AnnotateVCF, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.snp_filtered = "{0}/{1}.{2}.snp.filtered.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.final_vcf = "{0}/{1}.{2}.analysisReady.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.annotated_vcf = "{0}/{1}.{2}.analysisReady.annotated.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.annotated_vcf_gz = "{0}/{1}.{2}.analysisReady.annotated.vcf.gz".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.annotated_vcf_gz_index = "{0}/{1}.{2}.analysisReady.annotated.vcf.gz.tbi".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

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

        snpEff_cmd = ("{java} -Xmx5G -jar {snpEff} eff "
                     "GRCh37.74 -c {snpEff_cfg} "
                     "-interval {snpEff_interval} "
                     "-v -noMotif -noNextProt -noLog -nodownload -noStats "
                     "-o vcf {final_vcf}"
                     ).format(java=config().java,
                              snpEff=config().snpEff,
                              snpEff_cfg=config().snpEff_cfg,
                              snpEff_interval=config().snpEff_interval,
                              final_vcf=self.final_vcf)
        #In case of overwritting previous vcf.gz and vcf.gz.tbi files
        if os.path.isfile(self.annotated_vcf_gz):
            subprocess.check_call(['rm',self.annotated_vcf_gz])
        bgzip_cmd = ("{0} {1}").format(config().bgzip,self.annotated_vcf)
        tabix_cmd = ("{0} -f {1}").format(config().tabix,self.annotated_vcf_gz)

        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(snpEff_cmd + "\n")
                o.write(bgzip_cmd + "\n")
                o.write(tabix_cmd + "\n")
            with open(self.annotated_vcf,'w') as vcf_out, \
                open(self.log_file, "w") as log_fh:
                    p = subprocess.Popen(shlex.split(snpEff_cmd), stdout=vcf_out, stderr=log_fh)
                    p.wait()
            if p.returncode:
                raise subprocess.CalledProcessError(p.returncode, cmd)
            subprocess.check_call(shlex.split(snpEff_cmd))
            subprocess.check_call(shlex.split(bgzip_cmd))
            subprocess.check_call(shlex.split(tabix_cmd))

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="Dragen "+self.__class__.__name__,
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(CombineVariants)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class ArchiveSample(SGEJobTask):
    """Archive samples on Amplidata"""

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(ArchiveSample, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            config().scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.base_dir = "{0}/{1}/{2}.{3}".format(
            config().base,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.recal_bam = "{0}/{1}.{2}.realn.recal.bam".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.recal_bam_index = "{0}/{1}.{2}.realn.recal.bai".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.bam = "{0}/{1}.{2}.bam".format(
            self.base_dir,self.sample_name,self.pseudo_prepid)
        self.bam_index = "{0}/{1}.{2}.bam.bai".format(
            self.base_dir,self.sample_name,self.pseudo_prepid)
        self.annotated_vcf_gz = "{0}/{1}.{2}.analysisReady.annotated.vcf.gz".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.g_vcf_gz = "{0}/{1}.{2}.g.vcf.gz".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.g_vcf_gz_index = "{0}/{1}.{2}.g.vcf.gz.tbi".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.annotated_vcf_gz_index = "{0}/{1}.{2}.analysisReady.annotated.vcf.gz.tbi".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.copy_complete = "{0}/copy_complete".format(
            config().base)
        self.script_dir = "{0}/scripts".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)

        """
        self.cvg_binned = {0}/{1}.{2}.cvg_binned.format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.gq_binned = {0}/{1}.{2}.gq_binned.format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.align_metrics = {0}/{1}.{2}.alignment.metrics.txt.format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.cvg_metrics = {0}/{1}.{2}.cvg_metrics.txt.format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.cvg_metrics_ccds = {0}/{1}.{2}.cvg_metrics.ccds.txt.format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.cvg_metrics_X = {0}/{1}.{2}.cvg_metrics.X.txt.format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.cvg_metrics_Y = {0}/{1}.{2}.cvg_metrics.Y.txt.format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.dup_metrics = {0}/{1}.{2}.duplicates.txt.format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.variant_call_metrics = {0}/{1}.{2}.variant_calling_summary_metrics.format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.concordance_metrics = {0}/{1}.{2}.genotype_concordance_metrics.format(
        """

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

        """
        cmd = ("rsync -a --timeout=25000 -r "
              "{script_dir} {recal_bam} {recal_bam_index} {annotated_vcf_gz} "
              "{annotated_vcf_gz_index} {g_vcf_gz} {g_vcf_gz_index} "
              "{cvg_binned} {gq_binned} {align_metrics} 'cvg_metrics} "
              "{cvg_metrics_ccds} {cvg_metrics_X} {cvg_metrics_Y} {dup_metrics} "
              "{variant_call_metrics} {concordance_metrics} "
              "{base_dir}"
              ).format(self.__dict__))
        """

        cmd = ("rsync -a --timeout=25000 -r "
              "{script_dir} {recal_bam} {recal_bam_index} {annotated_vcf_gz} "
              "{annotated_vcf_gz_index} {g_vcf_gz} {g_vcf_gz_index} "
              "{base_dir}"
              ).format(recal_bam=self.recal_bam,
                      recal_bam_index=self.recal_bam_index,
                      script_dir=self.script_dir,
                      annotated_vcf_gz=self.annotated_vcf_gz,
                      annotated_vcf_gz_index=self.annotated_vcf_gz_index,
                      base_dir=self.base_dir,
                      g_vcf_gz=self.g_vcf_gz,
                      g_vcf_gz_index=self.g_vcf_gz_index,
                      sample_name=self.sample_name,
                      sample_type=self.sample_type.upper())

        try:
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)
            subprocess.check_call(['touch',self.copy_complete])
            """Original dragen BAM could be technically deleted earlier after the
            realigned BAM has been created on scratch space but it is safer to
            delete after the final realigned, recalculated BAM has been archived
            since our scratch space has failed in the past."""
            rm_cmd = ['rm',self.bam,self.bam_index]
            rm_folder_cmd = ['rm','-rf',self.scratch_dir]
            #subprocess.call(rm_cmd)
            #subprocess.call(rm_folder_cmd)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            cur.execute(INSERT_DRAGEN_STATUS.format(
                        sample_name=self.sample_name,
                        status="GATK Pipeline Complete",
                        DBID=DBID,
                        prepID=prepID,
                        pseudoPrepID=self.pseudo_prepid,
                        userID=userID))
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
        finally:
            if db.open:
                db.close()

    def requires(self):
        yield self.clone(AnnotateVCF)
        #yield self.clone(CvgBinning)
        #yield self.clone(GQBinning)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class SQLTarget(luigi.Target):
    """Target describing verification of the entries in the database"""
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


class MyExtTask(luigi.ExternalTask):
    """
    Checks whether the file specified exists on disk
    """

    file_loc = luigi.Parameter()
    def output(self):
        return luigi.LocalTarget(self.file_loc)


class CreateGenomeBed(SGEJobTask):
    """
    Use bedtools to create the input bed file for the coverage binning script
    """

    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter()

    ## System Parameters
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/rp2801/genome_cvg_temp/"

    def __init__(self,*args,**kwargs):
        super(CreateGenomeBed,self).__init__(*args,**kwargs)

        if not os.path.isdir(self.scratch): ## Recursively create the directory if it doesnt exist
            os.makedirs(self.scratch)

        self.sample_dir = os.path.join(self.scratch,self.sample_name)
        self.output_dir = os.path.join(self.sample_dir,'cvg_binned')
        self.genomecov_bed = os.path.join(self.sample_dir,self.sample_name+'.genomecvg.bed')
        self.recal_bam = os.path.join(self.sample_dir,self.sample_name+'.realn.recal.bam')
        self.genomecov_cmd = "{0} genomecov -bga -ibam {1} > {2}".format(
                              config().bedtools,self.recal_bam,
                              self.genomecov_bed)

        self.gzip_cmd = "bgzip {0}".format(self.genomecov_bed)
        self.gzipped_genomecov_bed = self.genomecov_bed+'.gz'
        self.tabix_cmd = "/nfs/goldstein/software/bin/tabix -p bed {0}".format(self.gzipped_genomecov_bed)
        self.tbi_file = self.gzipped_genomecov_bed+'.tbi'

        self.human_chromosomes = []
        self.human_chromosomes.extend(range(1, 23))
        self.human_chromosomes = [str(x) for x in self.human_chromosomes]
        self.human_chromosomes.extend(['X', 'Y','MT'])

    def requires(self):
        """
        The dependency is the presence of the bam file
        """

        return self.clone(PrintReads)

    def output(self):
        """
        Output is the genomecov bed file
        """

        return self.genomecov_bed

    def work(self):
        """
        Run the bedtools cmd
        """
        os.system(self.genomecov_cmd)
        #os.system(self.gzip_cmd)
        #os.system(self.tabix_cmd)

class CvgBinning(SGEJobTask):
    """
    Task to run the binning script for Coverage
    """

    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter()
        
    
    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/rp2801/genome_cvg_temp"
    
    
    def __init__(self,*args,**kwargs):

        super(CvgBinning,self).__init__(*args,**kwargs)
        
        self.sample_dir = os.path.join(self.scratch,self.sample_name)
        
        self.output_dir = os.path.join(self.sample_dir,'cvg_binned')
        self.genomecov_bed = os.path.join(self.sample_dir,self.sample_name+'.genomecvg.bed')
        
        if not os.path.isdir(self.output_dir): ## Recursively create the directory if it doesnt exist
            os.makedirs(self.output_dir)

        
        self.human_chromosomes = []
        self.human_chromosomes.extend(range(1, 23))
        self.human_chromosomes = [str(x) for x in self.human_chromosomes]
        self.human_chromosomes.extend(['X', 'Y','MT'])

        
        self.binning_cmd = "{0} {1} --sample_id {2} --block_size {3} --output_dir {4} {5} --mode coverage".format(
            config().pypy_loc,config().coverage_binner_loc,self.sample_name,config().block_size,self.output_dir,
            self.genomecov_bed)        

        
    def requires(self):
        """
        Dependency is the completion of the CreateGenomeBed Task
        """
                   
        return self.clone(CreateGenomeBed)
        
    
    def output(self):
        """
        Output from this task are 23 files for each chromosome
        """
        
        for chrom in self.human_chromosomes:
            file_loc = os.path.join(self.output_dir,self.sample_name+'_coverage_binned_'+config().block_size
                                +'_chr%s.txt'%chrom)
                
            yield luigi.LocalTarget(file_loc)
                

    def work(self):
        """ Run the binning script
        """

        os.system(self.binning_cmd)
        ## Delete the genomecvg bed file
        os.system('rm {0}'.format(self.genomecov_bed))

        
class GQBinning(SGEJobTask):
    """
    Task to run the binning script for GQ
    """

    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter()
    
    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/rp2801/genome_cvg_temp"
    
    
    def __init__(self,*args,**kwargs):

        super(GQBinning,self).__init__(*args,**kwargs)
       
        self.sample_dir = os.path.join(self.scratch,self.sample_name)
        self.output_dir = os.path.join(self.sample_dir,'gq_binned')
        self.gvcf = os.path.join(self.sample_dir,self.sample_name+'.g.vcf.gz') 

        if not os.path.isdir(self.output_dir): ## Recursively create the directory if it doesnt exist
            os.makedirs(self.output_dir)

        
        
        self.human_chromosomes = []
        self.human_chromosomes.extend(range(1, 23))
        self.human_chromosomes = [str(x) for x in self.human_chromosomes]
        self.human_chromosomes.extend(['X', 'Y','MT'])
                       
        self.binning_cmd = "{0} {1} --sample_id {2} --block_size {3} --output_dir {4} {5} --mode gq".format(
            config().pypy_loc,config().coverage_binner_loc,self.sample_name,config().block_size,self.scratch,
            self.gvcf)

        
    def requires(self):
        """
        Dependency is the completion of the CreateGenomeBed Task
        """            
            
        return self.clone(HaplotypeCaller)
        
    
    def output(self):
        """
        Output from this task are 23 files, one for each chromosome
        """
        
        for chrom in self.human_chromosomes:
            file_loc = os.path.join(self.scratch,self.sample_name+'_gq_binned_'+config().block_size
                                +'_chr%s.txt'%chrom)
                
            yield luigi.LocalTarget(file_loc)
                

    def work(self):
        """ Run the binning script
        """

        os.system(self.binning_cmd)

        
class AlignmentMetrics(SGEJobTask):
    """
    Run Picard AlignmentSummaryMetrics
    """

    base_directory = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    interval = luigi.Parameter()
    
    
    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"

    def __init__(self,*args,**kwargs):
        super(AlignmentMetrics,self).__init__(*args,**kwargs)
       
        self.sample_dir = os.path.join(self.scratch,self.sample_name)

        self.output_file_raw = os.path.join(self.sample_dir,self.sample_name+'.alignment.metrics.raw.txt')
        self.output_file = os.path.join(self.sample_dir,self.sample_name+'.alignment.metrics.txt')
        self.cmd = "{0} -XX:ParallelGCThreads=12 -jar {1} CollectAlignmentSummaryMetrics TMP_DIR={2} VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE={3} INPUT={4} OUTPUT={5}".format(config().java,config().picard,self.scratch,config().ref,self.bam,self.output_file_raw)               
        self.parser_cmd = """cat {0} | grep -v "^#" | awk -f ../sh/transpose.awk > {1}""".format(self.output_file_raw,self.output_file)
        
        
    def exists(self):
        """
        Check whether the output file is present
        """
        
        return self.clone(PrintReads)

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
        os.system(self.parser_cmd)

        
class CreateTargetFile(SGEJobTask):
    """ 
    Task for creating a new target file(the format used by Picard)
    from a bed file 
    """

    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/rp2801/genome_cvg_temp"

    def __init__(self,*args,**kwargs):
        super(CreateTargetFile,self).__init__(*args,**kwargs)
        
        self.bed_file = args[0]
        self.scratch = args[1]
        
        if not os.path.isdir(self.scratch): ## Recursively create the directory if it doesnt exist
            os.makedirs(self.scratch)

        self.bed_stem = utils.get_filename_from_fullpath(self.bed_file) + ".list"
        self.output_file = os.path.join(self.scratch,self.bed_stem)
        
        self.bedtointerval_cmd = ("{0} -jar {1} BedToIntervalList I={2} SD={3} OUTPUT={4}".format(config().java,config().picard,self.bed_file,config().seqdict_file,self.output_file))    

    
    
    def requires(self):
        """ 
        Dependency for this task is the existence of the bedfile
        and the human build37 sequence dictionary file
        """
        
        return [MyExtTask(config().bed_file),MyExtTask(config().seqdict_file)]


class RunCvgMetrics(SGEJobTask):
    """ 
    A luigi task
    """

    
    ## Task Specific Parameters  
    bam = luigi.Parameter()
    sample_name = luigi.Parameter()
    seq_type = luigi.Parameter()
    wgsinterval = luigi.BooleanParameter(default=False)
    bed_file = luigi.Parameter(default='/nfs/seqscratch11/rp2801/backup/coverage_statistics/ccds_r14.bed')
    x_y = luigi.BooleanParameter(default=True)
    create_targetfile = luigi.BooleanParameter(default=False)
    ## Full path to the temporary scratch directory
    ## (e.g./nfs/seqscratch11/rp2801/ALIGNMENT/exomes/als3445)
    scratch  = luigi.Parameter()
    
    ## System Parameters 
    n_cpu = 2
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/rp2801/genome_cvg_temp"

    def __init__(self,*args,**kwargs):
        """
        Initialize Class Variables
        :type args: List
        :type kwargs: Dict
        """

        super(RunCvgMetrics,self).__init__(*args,**kwargs)
           

        ## Define on the fly parameters
        self.sample_dir = os.path.join(self.scratch,self.sample_name)
        self.output_file = os.path.join(self.sample_dir,self.sample_name + ".cvg.metrics.txt")
        self.raw_output_file = os.path.join(self.sample_dir,self.sample_name + ".cvg.metrics.raw.txt")
        if self.x_y:
            self.output_file_X = os.path.join(self.sample_dir,self.sample_name+ ".cvg.metrics.X.txt")
            self.raw_output_file_X = os.path.join(self.sample_dir,self.sample_name + ".cvg.metrics.X.raw.txt")
            self.output_file_Y = os.path.join(self.sample_dir,self.sample_name+ ".cvg.metrics.Y.txt")
            self.raw_output_file_Y = os.path.join(self.sample_dir,self.sample_name + ".cvg.metrics.Y.raw.txt")

            
        if self.create_targetfile:
            self.bed_stem = utils.get_filename_from_fullpath(config().bed_file)
            self.target_file = os.path.join(self.sample_dir,self.bed_stem)

        else:
            self.target_file = config().target_file 

        
        ## Define shell commands to be run
        if self.seq_type.upper() == 'GENOME': ## If it is a genome sample
            if self.wgsinterval == True: ## Restrict wgs metrics to an interval
                self.cvg_cmd1 = """{0} -Xmx{1}g -XX:ParallelGCThreads=24 -jar {2} CollectWgsMetrics VALIDATION_STRINGENCY=LENIENT R={3} I={4} INTERVALS={5} O={6}MQ=20 Q=10""".format(config().java,config().max_mem,config().picard,config().ref,self.bam,self.target_file,self.raw_output_file)
                
            else: ## Run wgs metrics across the genome
                self.cvg_cmd1 = """{0} -Xmx{1}g -XX:ParallelGCThreads=24 -jar {2} CollectWgsMetrics VALIDATION_STRINGENCY=LENIENT R={3} I={4} O={5} MQ=20 Q=10""".format(config().java,config().max_mem,config().picard,config().ref,self.bam,self.raw_output_file)
               
                ## Run the cvg metrics on X and Y Chromosomes only
                self.cvg_sex1= """{0} -Xmx{1}g -XX:ParallelGCThreads=24 -jar {2} CollectWgsMetrics VALIDATION_STRINGENCY=LENIENT R={3} I={4} INTERVALS={5} O={6} MQ=20 Q=10"""
                self.cvg_cmd2 = self.cvg_sex1.format(config().java,config().max_mem,config().picard,config().ref,self.bam,config().bait_file_X,self.raw_output_file)              
                self.cvg_cmd3 = self.cvg_sex2.format(config().java,config().max_mem,config().picard,config().ref,self.bam,config().bait_file_Y,self.raw_output_file_X)              
                             
        else: ## Anything other than genome i.e. Exome or Custom Capture( maybe have an elif here to check further ?)
            self.cvg_cmd1 = """{0} -XX:ParallelGCThreads=24 -jar {1} CollectHsMetrics BI={2} TI={3} VALIDATION_STRINGENCY=LENIENT METRIC_ACCUMULATION_LEVEL=ALL_READS I={4} O={5} MQ=20 Q=10""".format(config().java,config().picard,config().bait_file,self.target_file,self.bam,self.raw_output_file)

            ## Run the Metrics cvg metrics on X and Y Chromosomes only
            self.cvg_sex2 =  """{0} -XX:ParallelGCThreads=24 -jar {1} CollectHsMetrics BI={2} TI={3} VALIDATION_STRINGENCY=LENIENT METRIC_ACCUMULATION_LEVEL=ALL_READS I={4} O={5} MQ=20 Q=10"""
            self.cvg_cmd2 =self.cvg_sex2.format(config().java,config().picard,config().bait_file,config().bait_file_X,self.bam,self.raw_output_file_X)
            self.cvg_cmd3 = self.cvg_sex2.format(config().java,config().picard,config().bait_file,config().bait_file_Y,self.bam,self.raw_output_file_Y)

        self.parser_cmd = """cat {0} | grep -v "^#" | awk -f /home/rp2801/git/dragen_pipe/dragen/python/transpose.awk > {1}"""
        self.parser_cmd1 = self.parser_cmd.format(self.raw_output_file,self.output_file)
        self.parser_cmd2 = self.parser_cmd.format(self.raw_output_file_X,self.output_file_X)
        self.parser_cmd3 = self.parser_cmd.format(self.raw_output_file_Y,self.output_file_Y)
        
       
            
    def output(self):
        """
        The output produced by this task 
        """
        yield luigi.LocalTarget(self.output_file)
        yield luigi.LocalTarget(self.output_file_X)
        yield luigi.LocalTarget(self.output_file_Y)
        

    def requires(self):
        """
        The dependency for this task is the CreateTargetFile task
        if the appropriate flag is specified by the user
        """

        if self.create_targetfile == True:
            yield [CreateTargetFile(self.bed_file,self.scratch),MyExtTask(self.bam)]
        else:
            yield MyExtTask(self.bam)
        
        
    def work(self):
        """
        Run Picard CalculateHsMetrics or WgsMetrics
        """
        ## Run the command
        os.system(self.cvg_cmd1)
        subprocess.check_output(self.parser_cmd1,shell=True)
        if self.x_y:
            os.system(self.cvg_cmd2)
            subprocess.check_output(self.parser_cmd2,shell=True)
            os.system(self.cvg_cmd3)
            subprocess.check_output(self.parser_cmd3,shell=True)

