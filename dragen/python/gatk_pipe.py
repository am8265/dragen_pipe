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
import re
import time
import glob
import warnings
import tabix

"""
Run samples through a luigized GATK pipeline after finishing the
Dragen based alignment
"""

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

def getUserID():
    """ Function to get useridi """

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
    """ Get location to the capturekit bedfile for a given
    pseudo_prepid

    pseudo_prepid : str : the pseudo_prepid value """
    db = get_connection("seqdb")
    try:
        cur = db.cursor()
        cur.execute(GET_CAPTURE_KIT_BED.format(
            pseudo_prepid=pseudo_prepid))
        capture_kit_bed = cur.fetchone()
    finally:
        if db.open:
            db.close()
    return capture_kit_bed

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

def getReadLength(pseudo_prepid):
    """ Get the chem version of all sequencing for pseudo_prepid and uses the
    most recent one.  There could be unlikely cases of samples sequenced
    with two different type of sequencing chemistry especially if we decide
    to do analysis with both exome and genome sequencing

    pseudo_prepid : str ; the pseudo_prepid value """

    db = get_connection("seqdb")
    try:
        cur = db.cursor()
        cur.execute(GET_READ_LENGTH.format(
            pseudo_prepid=pseudo_prepid))
        readLength = cur.fetchone()
    finally:
        if db.open:
            db.close()

    return readLength

class config(luigi.Config):
    """ Configuration class for instantiating parameters for this pipeline
    the values are read from luigi.cfg in the current folder """

    #programs
    bedtools = luigi.Parameter()
    bgzip = luigi.Parameter()
    gatk = luigi.Parameter()
    gatk4 = luigi.Parameter()
    java = luigi.Parameter()
    picard = luigi.Parameter()
    pypy = luigi.Parameter()
    snpEff = luigi.Parameter()
    snpSift = luigi.Parameter()
    tabix = luigi.Parameter()
    verifybamid = luigi.Parameter()
    vcffilter = luigi.Parameter()
    clineff = luigi.Parameter()
    #files
    ccds_bed_file = luigi.Parameter()
    coverage_binner = luigi.Parameter()
    gq_binner = luigi.Parameter()
    snpEff_cfg = luigi.Parameter()
    clineff_cfg = luigi.Parameter()
    target_file = luigi.Parameter()
    target_file_X = luigi.Parameter()
    target_file_Y = luigi.Parameter()
    transpose_awk = luigi.Parameter()
    subset_vcf_awk = luigi.Parameter()
    #gatk resources
    contam1000g_vcf = luigi.Parameter(description="The 1000g vcf file subsetted to SeqCap Regions to be used for contamination check")
    ref = luigi.Parameter()
    seqdict_file = luigi.Parameter()
    interval = luigi.Parameter()
    g1000 = luigi.Parameter()
    hapmap = luigi.Parameter()
    Mills1000g = luigi.Parameter()
    omni = luigi.Parameter()
    dbSNP = luigi.Parameter()
    annotatedbSNP = luigi.Parameter()
    exac = luigi.Parameter()
    #CNV Panel of Normals
    #roche_100_bp = luigi.Parameter()
    #roche_125_bp = luigi.Parameter()
    #genome_100_bp = luigi.Parameter()
    #genome_125_bp = luigi.Parameter()
    #variables
    create_targetfile = luigi.BooleanParameter()
    max_mem = luigi.IntParameter()
    block_size = luigi.Parameter(description='The block size over which to do the binning')
    cnv_target_padding = luigi.Parameter(description=("Num. of bp to pad each target during GATK's "
            "CalculateTargetCoverage CNV task"))
    #locations
    python_path = luigi.Parameter()
    #scratch = luigi.Parameter()
    base = luigi.Parameter()

class db(luigi.Config):
    """ Database config variable will be read from the
    db section of the config file """

    cnf = luigi.Parameter()
    seqdb_group = luigi.Parameter()
    dragen_group = luigi.Parameter()

class qcmetrics(luigi.Config):
    """ Parse in database field names for qc metrics
    from the qcmetrics section of the config file """

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

class CopyBam(SGEJobTask):
    """ Class for copying dragen aligned bam to a scratch location """

    sample_name = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    scratch = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(CopyBam, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        #self.base_dir = "{0}/{1}/{2}.{3}".format(
            #config().base,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)

        self.base_dir = "{0}/{1}/{2}.{3}".format("/nfs/fastq16/ALIGNMENT/BUILD37/DRAGEN",
                                                 self.sample_type.upper(),self.sample_name,
                                                 self.pseudo_prepid)
        self.base_log = "{0}/logs".format(self.base_dir)
        self.bam = "{0}/{1}.{2}.bam".format(self.base_dir,
                                            self.sample_name,self.pseudo_prepid)
        self.bam_index = "{0}/{1}.{2}.bam.bai".format(self.base_dir,
                                                      self.sample_name,self.pseudo_prepid)
        self.scratch_bam_index = os.path.join(self.scratch_dir,'{0}.{1}.bam.bai'.format(self.sample_name,self.pseudo_prepid))
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
        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
            db.close()

            os.system("mkdir -p {0}/scripts".format(self.scratch_dir))
            os.system("mkdir -p {0}/logs".format(self.scratch_dir))
            cmd = ("rsync -a --timeout=20000 -r {bam} {bam_index} {base_log} {scratch_dir}/"
                  ).format(bam=self.bam,
                           bam_index=self.bam_index,
                           base_log=self.base_log,
                           scratch_dir=self.scratch_dir)
            with open(self.script,'w') as o:
                o.write(cmd + "\n")

            DEVNULL = open(os.devnull, 'w')
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()
        finally:
            if db.open:
                db.close()

    def output(self):
        ## Need to have a convoluted logic here since there are samples which have their dragen alignment written
        ## to the scratch dir itself, so for samples which do not have a bam index in the scratch directory we will
        ## go ahead with the rsync and regular db update for copy bam

        ## For samples which have do have a dragen bam in the scratch directory, the work function will not be called
        ## and hence the copy bam task will be absent from the dragen_pipeline_step table, this needs to be fixed

        ## Note, for the future it is better if the dragen bam could be created in this step ,
        ## luigi tasks are supposed to be idempotant (http://datapipelinearchitect.com/luigi-force-rerun/)
        ## Since the dragen bam is deleted after the sample is archived the logical way to re-reun a pipeline
        ## would be deleting database table entries, but since the dragen bam can never be recovered again the
        ## pipeline would fail to execute

        if os.path.isfile(self.scratch_bam_index):
            return luigi.LocalTarget(self.scratch_bam_index)
        else:
            return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                             pipeline_step_id=self.pipeline_step_id)


class RealignerTargetCreator(SGEJobTask):
    """ Class for creating targets for indel realignment BAMs from Dragen """

    sample_name = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    scratch = luigi.Parameter()

    n_cpu = 4
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(RealignerTargetCreator, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
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
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
            db.close()
            
            with open(self.script,'w') as o:
                o.write(cmd + "\n")

            with open(os.devnull,'w') as DEVNULL:
                subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            
            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()
        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(CopyBam)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)


class IndelRealigner(SGEJobTask):
    """ class to create BAM with realigned BAMs
    """
    
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(IndelRealigner, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
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

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
        
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
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
            db.close()
            with open(self.script,'w') as o:
                o.write(cmd + "\n")
                
            with open(os.devnull,'w') as DEVNULL:
                subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()

        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(RealignerTargetCreator)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class BaseRecalibrator(SGEJobTask):
    """ Class to create a recalibration table with realigned BAMs
    """
    
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()

    n_cpu = 4
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(BaseRecalibrator, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.realn_bam = "{0}/{1}.{2}.realn.bam".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.recal_table = "{0}/{1}.{2}.recal_table".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
        
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T BaseRecalibrator "
            "-L {interval} "
            "-I {realn_bam} "
            "-o {recal_table} "
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
                dbSNP=config().dbSNP)

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
            db.close()
            with open(self.script,'w') as o:
                o.write(cmd + "\n")
                
            with open(os.devnull, 'w') as DEVNULL:
                subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()

        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(IndelRealigner)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class PrintReads(SGEJobTask):
    """ Class to create BAM with realigned and recalculated BAMs
    """
    
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    
    n_cpu = 4
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(PrintReads, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
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

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
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
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
            db.close()
            
            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            with open(os.devnull, 'w') as DEVNULL:
                subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()

            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()
        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(BaseRecalibrator)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)


class HaplotypeCaller(SGEJobTask):
    """ Class to create gVCFs
    """
    
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()

    n_cpu = 4
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(HaplotypeCaller, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
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

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
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
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
            db.close()
            
            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            with open(os.devnull, 'w') as DEVNULL:
                subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            rm_cmd = "rm {0} {1}".format(self.realn_bam,self.scratch_bam)
            if os.path.exists(self.realn_bam) and os.path.exists(self.scratch_bam):
                with open(self.script,'a') as o:
                    o.write(rm_cmd + "\n")
                with open(self.log_file, "a") as log_fh, open(os.devnull,'w') as DEVNULL:
                    p = subprocess.Popen(shlex.split(rm_cmd), stdout=DEVNULL, stderr=log_fh)
                    p.wait()                    
                    if p.returncode:
                        raise subprocess.CalledProcessError(p.returncode, rm_cmd)
            
            
            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()

            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()
        finally:
            if db.open:
                db.close()
                
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
    scratch = luigi.Parameter()
    
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(GenotypeGVCFs, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
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

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
            db.close()
        except MySQLdb.Error, e:
            raise Exception("ERROR %d IN CONNECTION: %s" % (e.args[0], e.args[1]))
        finally:
            if db.open:
                db.close()

    def work(self):
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
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
            db.close()
            
            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            with open(os.devnull, 'w') as DEVNULL:
                subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()

            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()
            
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
    scratch = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(SelectVariantsSNP, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
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

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
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
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
            db.close()

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            with open(os.devnull, 'w') as DEVNULL:
                subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()

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
    scratch = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(SelectVariantsINDEL, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
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
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
            db.close()

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            with open(os.devnull, 'w') as DEVNULL:
                subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()

        finally:
            if db.open:
                db.close()

    def requires(self):
        if self.sample_type.lower() == 'exome':
            return self.clone(ApplyRecalibrationSNP)
        if self.sample_type.lower() == 'genome':
            return self.clone(ApplyRecalibrationSNP)
        elif self.sample_type.lower() == 'custom_capture':
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
    scratch = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(VariantRecalibratorSNP, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
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

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):

        ## Running both Exomes and Genomes the same way, i.e. we will exlcude DP, omni is set to be a truth set, also added in ExAC SNPs as a training
        ## set which can contain both TP and FP with the same prior likelihood as the 1000G training set. 
        ## See these links : 
        ## For exomes : 
        ## https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php
        ## For genomes : 
        ## https://software.broadinstitute.org/gatk/guide/article?id=1259

        cmd = ("{java} -Xmx{max_mem}g "
               "-jar {gatk} "
               "-R {ref} "
               "-T VariantRecalibrator "
               "-L {interval} "
               "--log_to_file {log_file} "
               "--input {snp_vcf} "
               "-an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum "
               "-mode SNP "
               "--maxGaussians 4 "
               "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 "
               "-recalFile {snp_recal} "
               "-tranchesFile {snp_tranches} "
               "-rscriptFile {snp_rscript} "
               "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} "
               "-resource:omni,known=false,training=true,truth=true,prior=12.0 {omni} "
               "-resource:1000G,known=false,training=true,truth=false,prior=10.0 {g1000} "
               "-resource:ExAc,known=false,training=true,truth=false,prior=10.0 {exac} "
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
                            g1000=config().g1000,
                            exac=config().exac)

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
            db.close()

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            with open(os.devnull, 'w') as DEVNULL:
                subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()

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
    scratch = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(VariantRecalibratorINDEL, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
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

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
        cmd = ("{java} -Xmx{max_mem}g "
            "-jar {gatk} "
            "-R {ref} "
            "-T VariantRecalibrator "
            "-L {interval} "
            "--log_to_file {log_file} "
            "--input {indel_vcf} "
            "-an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum "
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
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
            db.close()

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            with open(os.devnull, 'w') as DEVNULL:
                subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()

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
    scratch = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(ApplyRecalibrationSNP, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
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

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
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
               "--ts_filter_level 90.0 "  ## We will be using a sensitvity filter of 90.0 to get the 90.0to99.0 tranche as well. The resulting vcf will have 4 stratifications , PASS , CONFIDENT, INTERMEDIATE, NOT CONFIDENT
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
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
            db.close()
            
            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            with open(os.devnull, 'w') as DEVNULL:
                subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()
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
    scratch = luigi.Parameter()
    
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(ApplyRecalibrationINDEL, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
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

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
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
            "--ts_filter_level 90.0 "
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
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            db.commit()
            db.close()
            
            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            with open(os.devnull, 'w') as DEVNULL:
                subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()
            
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
    scratch = luigi.Parameter()
    
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(VariantFiltrationSNP, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.snp_vcf = "{0}/{1}.{2}.snp.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.snp_filtered = "{0}/{1}.{2}.snp.filtered.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
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
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            db.commit()
            db.close()
            
            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            with open(os.devnull, 'w') as DEVNULL:
                subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()

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
    scratch = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(VariantFiltrationINDEL, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.indel_vcf = "{0}/{1}.{2}.indel.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.indel_filtered = "{0}/{1}.{2}.indel.filtered.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
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
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            db.commit()
            db.close()

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            with open(os.devnull, 'w') as DEVNULL:
                subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)

            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()

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
    scratch = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(CombineVariants, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
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
        self.final_vcf = "{0}/{1}.{2}.analysisReady.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        
        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
        """Merges SNP and INDEL vcfs.  Using the filtered SNP vcf header as a
        base, variant type/sample type specific ##FILTERs are added to the header.
        Additionally the AnnoDBID annotation is added.  After finishing reading
        through the header the SNP vcf is read again with only variants
        being outputted.  The same happens with the INDEL vcf.  The resulting
        vcf is then sorted producing the analysisReady.vcf.
        """
        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            db.commit()
            db.close()
            
            filter_flag = 0
            info_flag = 0
            
            with open(self.tmp_vcf,'w') as vcf_out, open(self.snp_filtered) as header:
                for line in header.readlines():
                    if line[0] == '#':
                        if line[0:8] == '##FILTER' and filter_flag == 0 :
                            filter_flag = 1
                            #SNP specific FILTER whether using VQSR or snp filtering
                            if self.sample_type == 'exome'.lower() or self.sample_type.lower() =='genome':            ## Will a 90.0to99.00 tranche be needed in the header ?
                                vcf_out.write('##FILTER=<ID=VQSRTrancheSNP90.00to99.00,Description="Truth sensitivity tranche level for SNP model"\n')
                                vcf_out.write('##FILTER=<ID=VQSRTrancheSNP99.00to99.90,Description="Truth sensitivity tranche level for SNP model"\n')
                                vcf_out.write('##FILTER=<ID=VQSRTrancheSNP99.90to100.00+,Description="Truth sensitivity tranche level for SNP model"\n')
                                vcf_out.write('##FILTER=<ID=VQSRTrancheSNP99.90to100.00,Description="Truth sensitivity tranche level for SNP model"\n')
                            else:
                                vcf_out.write('##FILTER=<ID=SNP_filter,Description="QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0">\n')

                            #Indel specific filters whether using VQSR or indel filtering
                            if self.sample_type.lower() =='genome':
                                vcf_out.write('##FILTER=<ID=VQSRTrancheINDEL90.00to99.00,Description="Truth sensitivity tranche level for INDEL model"\n')
                                vcf_out.write('##FILTER=<ID=VQSRTrancheINDEL99.00to99.90,Description="Truth sensitivity tranche level for INDEL model"\n')
                                vcf_out.write('##FILTER=<ID=VQSRTrancheINDEL99.90to100.00+,Description="Truth sensitivity tranche level for INDEL model"\n')
                                vcf_out.write('##FILTER=<ID=VQSRTrancheINDEL99.90to100.00,Description="Truth sensitivity tranche level for INDEL model"\n')
                            else:
                                vcf_out.write('##FILTER=<ID=INDEL_filter,Description="QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0">\n')

                        #AnnoDBID annotation will be added during the annotation pipeline
                        if line[0:13] == "##INFO=<ID=AN":
                                vcf_out.write('##INFO=<ID=AnnoDBID,Number=1,Type=String,Description="AnnoDBID">\n')
                                
                        if line[0:8] == '##INFO' and info_flag == 0:
                            info_flag = 1
                            ## Adding this wierd SAO info tag since there are rare cases where a variant has this information, but it is absent in the header
                            ## causing read backed phasing to crash.
                            vcf_out.write("""##INFO=<ID=SAO,Number=1,Type=Integer,Description="Variant Allele Origin: 0. unspecified, 1. Germline, 2. Somatic, 3. Both">\n""")
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

            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()

        finally:
            if db.open:
                db.close()

    def requires(self):
        if self.sample_type.lower() == 'exome':
            return self.clone(VariantFiltrationINDEL)
        elif self.sample_type.lower() == 'custom_capture':
            return self.clone(VariantFiltrationINDEL)
        elif self.sample_type.lower() == 'genome':
            return self.clone(ApplyRecalibrationINDEL)
        else:
            raise Exception, "Sample type: %s not supported in this module" % self.sample_type

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class RBP(SGEJobTask):
    """ Run GATK's Read Backed Phasing"""

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(RBP,self).__init__(*args,**kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.recal_bam = "{0}/{1}.{2}.realn.recal.bam".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.deannotated_vcf = "{0}/{1}.{2}.analysisReady.deannotated.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.deannotated_vcf_gz = "{0}/{1}.{2}.analysisReady.deannotated.vcf.gz".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.annotated_vcf_gz = "{0}/{1}.{2}.analysisReady.annotated.vcf.gz".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.final_vcf = "{0}/{1}.{2}.analysisReady.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.phased_vcf = "{0}/{1}.{2}.analysisReady.phased.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.phased_vcf_gz = "{0}/{1}.{2}.analysisReady.phased.vcf.gz".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                    step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
        """ Run the actual command """

        self.bgzip_cmd = "{0} -f {1}".format(config().bgzip,self.final_vcf)
        self.final_vcf_gz = self.final_vcf+".gz"
        self.tabix_cmd = "{0} -f {1}".format(config().tabix,self.final_vcf_gz)
        if os.path.exists(self.final_vcf_gz):
            self.final_vcf = self.final_vcf_gz

        self.cmd = (
            """ {java} -jar {gatk} -T ReadBackedPhasing -R {ref} """
            """ -I {recal_bam} --variant {vcf} -o {phased_vcf} """
            """ --phaseQualityThresh 20.0 --maxGenomicDistanceForMNP 2 """
            """ --enableMergePhasedSegregatingPolymorphismsToMNP -U ALLOW_SEQ_DICT_INCOMPATIBILITY >& {log}""".format(
                java=config().java,gatk=config().gatk,ref=config().ref,
                recal_bam=self.recal_bam,vcf=self.final_vcf,
                phased_vcf=self.phased_vcf,log=self.log_file))

        with open(self.script,'w') as OUT:
                OUT.write(self.cmd)

        run_shellcmd("seqdb",self.pipeline_step_id,self.pseudo_prepid,
                     [self.cmd],self.sample_name,self.__class__.__name__)

    def requires(self):
        """ The task dependency """

        return self.clone(CombineVariants)

    def output(self):
        """ The output from this task """

        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

class FixMergedMNPInfo(SGEJobTask):
    """ The MNPs phased by RBP are missing the DP,AD and other annotation info
    this task which check variants which are missing these information, fetch
    the value for corresponding to that variant site from the vcf prior to RBP
    and add these information to the fixed vcf """

    sample_name = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self,*args,**kwargs):
        super(FixMergedMNPInfo,self).__init__(*args,**kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
        self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
        self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.phased_vcf = "{0}/{1}.{2}.analysisReady.phased.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.fixed_vcf = "{0}/{1}.{2}.analysisReady.fixed.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.final_vcf = "{0}/{1}.{2}.analysisReady.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

    def requires(self):
        return self.clone(RBP)

    def work(self):
        self.LOG_FILE = open(self.log_file,'w')
        tries = 5
        for i in range(tries):
            try:
                db = get_connection("seqdb")
                cur = db.cursor()
                cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                    pseudo_prepid=self.pseudo_prepid,
                    pipeline_step_id=self.pipeline_step_id))
                db.commit()
                db.close()

                ## bgzip and tabix index the final vcf file
                self.bgzip_cmd = "{0} -f {1}".format(config().bgzip,self.final_vcf)
                self.final_vcf_gz = self.final_vcf+".gz"
                self.tabix_cmd = "{0} -f {1}".format(config().tabix,self.final_vcf_gz)
                if os.path.exists(self.final_vcf):
                    os.system(self.bgzip_cmd)
                os.system(self.tabix_cmd)
                self.fix_phased_vcf(self.phased_vcf,self.final_vcf_gz,self.fixed_vcf)
                self.LOG_FILE.close()

                ## Remove temp files
                if os.path.isfile(self.phased_vcf):
                    subprocess.check_call(['rm',self.phased_vcf])

                db = get_connection("seqdb")
                cur = db.cursor()
                update_dragen_status(self.pseudo_prepid,self.sample_name,
                                     self.__class__.__name__,cur)
                cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                    pseudo_prepid=self.pseudo_prepid,
                    pipeline_step_id=self.pipeline_step_id))
                db.commit()
                db.close()

            except MySQLdb.Error, e:
                if i == tries - 1:
                    raise Exception("ERROR %d IN CONNECTION: %s"%(e,args[0],e.args[1]))
                else:
                    time.sleep(60)
                    continue
            finally:
                if db.open:
                    db.close()
                self.LOG_FILE.close()
                break

    def fix_phased_vcf(self,phased_vcf,deannotated_vcf_gz,fixed_vcf):
        """ Fix the missing DP and AD fields for the phased variants

        phased_vcf : str ; path to the phased vcf file
        deannotated_vcf : str ; path to the deannotated vcf file(should
        be gzipped and tabix indexed
        """

        with open(phased_vcf,'r') as PHASED_VCF, open(fixed_vcf,'w') as OUT:
            for line in PHASED_VCF:
                line = line.strip('\n')
                if line[0] == '#':
                    print >> OUT,line
                elif self.check_dp_info(line): ## Check if DP is present in the info field
                    print >> OUT,line
                else:
                    ## Check for a MNP
                    if self.check_mnp(line):
                        ## get DP,AD_ref,AD_alt
                        chrom = line.split('\t')[0]
                        pos = int(line.split('\t')[1])
                        ad,dp = self.get_ad_dp_from_vcf(deannotated_vcf_gz,chrom,pos)
                        annotations = self.get_gatk_annotations_from_vcf(
                            deannotated_vcf_gz,chrom,pos)
                        ## Construct a new vcf line
                        new_line = self.construct_fixed_vcf_line(line,[ad,dp,annotations])
                        print new_line
                        print >> OUT,new_line
                    else: ## A non MNP with no DP, raise an exception
                        raise Exception("A non MNP encountered with no DP in info field")

    def construct_fixed_vcf_line(self,line,missing_vals):
        """ Construct the fixed line 

        line : str ; the vcf line
        missing_vals : [str,str,str] ; the missing dp,ad and
        gatk annotations like ExcessHet,FS,VQSLOD,etc.

        returns : the line with the missing values fixed
        """

        info_fields = line.split('\t')[-2].split(':')
        vals = line.split('\t')[-1].split(':')

        ## Insert the new tags at the correct position
        info_fields.insert(1,'AD')
        info_fields.insert(2,'DP')

        vals.insert(1,missing_vals[0])
        vals.insert(2,missing_vals[1])

        new_info = ':'.join(info_fields)
        new_vals = ':'.join(vals)

        ## Fix the missing GATK annotations
        annotations = line.split('\t')[7]
        new_annotations = annotations + ';' + missing_vals[2]

        contents = line.split('\t')
        contents[-1] = new_vals
        contents[-2] = new_info
        contents[7] = new_annotations

        ## Return the fixed line
        return '\t'.join(contents)

    def check_dp_info(self,line):
        """ Checks whether dp is present in the info field
        returns True if DP is present, otherwise returns False

        line : str ; a non header line from the vcf file

        returns : Bool 
        """

        return ("DP" in line.split('\t')[-2].split(':'))

    def check_mnp(self,line):
        """ Checks whether the variant is a MNP
        returns True if it is, otherwise returns False

        line : str ; a non header line from the vcf file

        returns : Bool 
        """

        def check_equality(r,a):
            return (len(r) == len(a))

        ref = line.split('\t')[3]
        alt = line.split('\t')[4]
        if len(ref) == 1:
            return False
        else:
            if len(alt.split(',')) > 1: ## A multiallelic site
                a1,a2 = alt.split(',')
                if check_equality(ref,a1) or check_equality(ref,a2):
                    raise Exception("A MNP at a multi allelic site !")
            else:
                return check_equality(ref,alt)                    

    def get_gatk_annotations_from_vcf(self,vcf,chrom,pos):
        """ Get GATK specific annotations like FS,VQSLOD,ExcessHet
        which are absent in MNPs after RBP has been run

        vcf : str ; path to the vcf file
        chrom : str ; the chromosome number
        pos : int ; the position for the variant

        returns : The annotations with the AC,AF,AN removed
        """

        try:
            tb = tabix.open(vcf)
        except IOError as (errno,etrerror):
            print >> self.LOG_FILE,"Error with tabix, could not open the vcf file"
            " check to see if it is indexed properly"
        except:
            print >> self.LOG_FILE,"Something wrong with opening the indexed vcf file"
            sys.exit(1)

        return (';'.join(list(tb.query(chrom,pos,pos))[0][7].split(';')[3:]))

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
        print >> self.LOG_FILE,message
        raise Exception(message)

    def get_ad_dp_from_vcf(self,vcf,chrom,pos):
        """ Get DP, AD_REF, AD_ALT
        from a VCF file

        vcf : str ; path to the vcf file
        chrom : str ; the chromosome number
        pos : int ; the position for the variant

        returns : [dp,ad_ref,ad_alt]
        """

        try:
            tb = tabix.open(vcf)
        except IOError as (errno,strerror):
            print >> self.LOG_FILE,"Error with tabix, could not open the vcf file"
            " check to see if it is indexed properly"
            sys.exit(1)
        except:
            print >> self.LOG_FILE,"Something wrong with opening the indexed vcf file"
            sys.exit(1)

        records = list(tb.query(chrom,pos,pos))
        ## The above query returns all variant sites
        ## spanning the pos which we have queried for
        ## We only want the exact pos site, so check
        ## for that condition
        record = self.get_correct_record(records,pos)
        info_fields = record[-2].split(':')
        values = record[-1].split(':')
        found = 0
        for i in range(0,len(info_fields)):
            if info_fields[i] == 'DP':
                dp = values[i]
                found+=1
            if info_fields[i] == 'AD':
                ad = values[i]  ## How to deal with multi allelic sites ? 
                found+=1

        if found != 2:
            raise Exception("Could not parse the vcf to get DP/AD info for "
                            "the phased variants")

        return [ad,dp]

class AnnotateVCF(SGEJobTask):
    """ Annotate the final vcf with SNPEff
    """

    sample_name = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(AnnotateVCF, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.{3}.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        self.snp_filtered = "{0}/{1}.{2}.snp.filtered.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.fixed_vcf = "{0}/{1}.{2}.analysisReady.fixed.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.annotated_vcf = "{0}/{1}.{2}.analysisReady.annotated.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.annotated_vcf_gz = "{0}/{1}.{2}.analysisReady.annotated.vcf.gz".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.annotated_vcf_gz_index = "{0}/{1}.{2}.analysisReady.annotated.vcf.gz.tbi".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):


        annotation_cmd = (""" /nfs/goldstein/software/jdk1.8.0_05/bin/java """
                          """ -Xmx6g -jar {0} -c {1} -v -db """
                          """ {2}"""
                          """ GRCh37.87 {3} > {4} 2> {5} """.format(
                              config().clineff,config().clineff_cfg,
                              config().annotatedbSNP,self.fixed_vcf,
                              self.annotated_vcf,self.log_file))

        bgzip_cmd = ("{0} -f {1}").format(config().bgzip,self.annotated_vcf)
        tabix_cmd = ("{0} -f {1}").format(config().tabix,self.annotated_vcf_gz)

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))

            db.commit()
            db.close()

            with open(self.script,'w') as o:
                o.write(annotation_cmd + "\n")
                o.write(bgzip_cmd + "\n")
                o.write(tabix_cmd + "\n")

            proc = subprocess.Popen(annotation_cmd,shell=True)
            proc.wait()
            if proc.returncode:
                raise subprocess.CalledProcessError(proc.returncode,annotation_cmd)
            subprocess.check_call(shlex.split(bgzip_cmd))
            subprocess.check_call(shlex.split(tabix_cmd))


            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()

            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()
        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(FixMergedMNPInfo)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)


class CleanDirectory(SGEJobTask):
    """ For samples which have already been archived , the pipeline will be rerun on fastq16 itself,
    in that case run this task before Archive to clear the temporary files from the directory """

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()

    ## System Parameters
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self,*args,**kwargs):
        super(CleanDirectory,self).__init__(*args,**kwargs)
        ## The scratch and base will be the same directories if this task
        ## is running
        self.seqtype_dir = os.path.join(self.scratch,self.sample_type.upper())
        self.sample_dir = os.path.join(self.seqtype_dir,self.sample_name+'.'+self.pseudo_prepid)
        kwargs["sample_name"] = self.sample_name
        super(CleanDirectory, self).__init__(*args,**kwargs)

        ## The intermediate files created in the pipeline run (note: many of these are wild cards)
        self.sample_stem = self.sample_name+'.'+self.pseudo_prepid
        self.analysis_ready_vcf = os.path.join(self.sample_dir,self.sample_stem+'.analysisReady.vcf*')
        self.indel_vcf = os.path.join(self.sample_dir,self.sample_stem+'.indel*vcf*')
        self.raw_vcf = os.path.join(self.sample_dir,self.sample_stem+'.raw.vcf*')
        self.snp_files = os.path.join(self.sample_dir,self.sample_stem+'*snp*')
        self.tmp_vcf = os.path.join(self.sample_dir,self.sample_stem+'.tmp.vcf')
        self.eval_vcf = os.path.join(self.sample_dir,'eval.vcf*')
        self.truth_vcf = os.path.join(self.sample_dir,'truth.vcf*')
        self.temp_geno_file = os.path.join(self.sample_dir,'temp_geno_concordance.txt')
        self.targets = os.path.join(self.sample_dir,'targets.interval')
        self.raw_cvg_files = os.path.join(self.sample_dir,self.sample_stem,'*raw*')
        self.genomecvg = os.path.join(self.sample_dir,self.sample_name+'.genomecvg.bed')
        self.dragen_bam = os.path.join(self.sample_dir,self.sample_stem+'.bam*')
        self.recal_table = os.path.join(self.sample_dir,self.sample_stem+'.recal_table')

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
        cmd = ("rm {analysis_ready_vcf} {indel_vcf} {raw_vcf} {snp_files} {tmp_vcf} {eval_vcf} "
               "{truth_vcf} {temp_geno_file} {targets} {raw_cvg_files} {genomecvg} {dragen_bam} {recal_table}"
              ).format(**self.__dict__)

        run_shellcmd("seqdb",self.pipeline_step_id,self.pseudo_prepid,
                     [cmd],self.sample_name,self.__class__.__name__)

    def requires(self):
        """ The requirement for this task is the presence of the analysis ready
        vcf file from the GATK pipeline """
        yield self.clone(AnnotateVCF)
        yield self.clone(GQBinning)
        yield self.clone(CvgBinning)
        yield self.clone(UpdateSeqdbMetrics)

    def output(self):
        """ The result from this task is the creation of the metrics file """
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

class SubsetVCF(SGEJobTask):
    """ For SRR/Any future deemed to be low quality sequenced samples, subset the vcf to only the capturekit regions """

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    ## System Parameters
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self,*args,**kwargs):
        super(SubsetVCF,self).__init__(*args,**kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
        self.annotated_vcf_gz = "{0}/{1}.{2}.analysisReady.annotated.vcf.gz".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.original_vcf_gz = "{0}/{1}.{2}.analysisReady.annotated.original.vcf.gz".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.restricted_vcf = "{0}/{1}.{2}.analysisReady.annotated.restricted.vcf".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.restricted_vcf_gz = "{0}/{1}.{2}.analysisReady.annotated.restricted.vcf.gz".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.log_file = "{0}/logs/{1}.{2}.subsetvcf.log".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.subset_cmd = "{0} {1} -B {2} -h > {3} 2> {4}".format(
            config().tabix,self.annotated_vcf_gz,self.capture_kit_bed,
            self.restricted_vcf,self.log_file)
        self.bgzip_cmd = "{0} -f {1}".format(config().bgzip,self.restricted_vcf)
        self.rename1 = "mv {0} {1}".format(self.annotated_vcf_gz,self.original_vcf_gz)
        self.rename2 = "mv {0} {1}".format(self.restricted_vcf_gz,self.annotated_vcf_gz)
        self.tabix_cmd1 = "{0} -f {1}".format(config().tabix,self.original_vcf_gz)
        self.tabix_cmd2 = "{0} -f {1}".format(config().tabix,self.annotated_vcf_gz)
        self.script = "{0}/scripts/{1}.{2}.{3}.sh".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)

        with open(self.script,'w') as OUT:
            OUT.write(self.subset_cmd+'\n')
            OUT.write(self.bgzip_cmd+'\n')
            OUT.write(self.rename1+'\n')
            OUT.write(self.rename2+'\n')
            OUT.write(self.tabix_cmd1+'\n')
            OUT.write(self.tabix_cmd2+'\n')

        ## Get the pipeline step id
        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
        run_shellcmd("seqdb",self.pipeline_step_id,self.pseudo_prepid,
                     [self.subset_cmd,self.bgzip_cmd,self.rename1,self.rename2,self.tabix_cmd1,self.tabix_cmd2,],self.sample_name,
                     self.__class__.__name__)

    def requires(self):
        return self.clone(AnnotateVCF)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

class ArchiveSample(SGEJobTask):
    """ Archive samples on Amplidata """

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self, *args, **kwargs):
        super(ArchiveSample, self).__init__(*args, **kwargs)
        self.scratch_dir = "{0}/{1}/{2}.{3}".format(
            self.scratch,self.sample_type.upper(),self.sample_name,self.pseudo_prepid)
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
        self.original_vcf_gz = "{0}/{1}.{2}.analysisReady.annotated.original.vcf.gz".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.original_vcf_gz_index = "{0}/{1}.{2}.analysisReady.annotated.original.vcf.gz.tbi".format(
            self.scratch_dir,self.sample_name,self.pseudo_prepid)
        self.copy_complete = "{0}/copy_complete".format(
            config().base)
        self.script_dir = "{0}/scripts".format(
            self.scratch_dir)
        self.log_dir = "{0}/logs".format(
            self.scratch_dir)

        ## qc metric files 
        self.sample_dir = self.scratch_dir
        self.alignment_out = os.path.join(self.sample_dir,"{0}.{1}.alignment.metrics.txt".format(self.sample_name,self.pseudo_prepid))
        self.cvg_out = os.path.join(self.sample_dir,"{0}.{1}.cvg.metrics.txt".format(self.sample_name,self.pseudo_prepid))
        self.cvg_ccds_out = os.path.join(self.sample_dir,"{0}.{1}.cvg.metrics.ccds.txt".format(self.sample_name,self.pseudo_prepid))
        self.cvg_X_out = os.path.join(self.sample_dir,"{0}.{1}.cvg.metrics.X.txt".format(self.sample_name,self.pseudo_prepid))
        self.cvg_Y_out = os.path.join(self.sample_dir,"{0}.{1}.cvg.metrics.Y.txt".format(self.sample_name,self.pseudo_prepid))
        if self.sample_type.upper() != 'GENOME': ## exome and custom captures will have capture specificity as well
            self.cvg_cs_out = os.path.join(self.sample_dir,"{0}.{1}.cvg.metrics.cs.txt".format(self.sample_name,self.pseudo_prepid))
        self.dup_out = os.path.join(self.sample_dir,"{0}.{1}.duplicates.txt".format(self.sample_name,self.pseudo_prepid))
        self.variant_call_out = os.path.join(self.sample_dir,"{0}.{1}.variant_calling_summary_metrics".format(self.sample_name,self.pseudo_prepid))
        self.geno_concordance_out = os.path.join(self.sample_dir,"{0}.{1}.genotype_concordance_metrics".format(self.sample_name,self.pseudo_prepid))
        self.contamination_out = os.path.join(self.sample_dir,"{0}.{1}.contamination.selfSM".format(self.sample_name,self.pseudo_prepid))
        self.cvg_binned = os.path.join(self.sample_dir,'cvg_binned')
        self.gq_binned = os.path.join(self.sample_dir,'gq_binned')

        ## The base directory will need to be modified if it is a pgm sample
        self.pgm_base = '/nfs/pgmclin/ALIGNMENT/BUILD37/DRAGEN/'
        if (self.sample_name.upper().startswith('PGMCLIN') or
            self.sample_name.upper().startswith('PGMVIP')):
            self.base_dir = "{0}/{1}/{2}.{3}".format(
                self.pgm_base,self.sample_type.upper(),self.sample_name,
                self.pseudo_prepid)

        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=self.__class__.__name__))
            self.pipeline_step_id = cur.fetchone()[0]
        finally:
            if db.open:
                db.close()

    def work(self):
        ## copy over the fastq directory as well
        self.fastq_loc = os.path.join(self.sample_dir,'fastq')
        ## create a dummy one if it doesnt exist
        if not os.path.exists(self.fastq_loc):
            os.makedirs(self.fastq_loc)
        ## create a dummy geno_concordance file if it doesnt exist , this is because there are times when this task
        ## will not be scheduled (failure to locate a production vcf) , then the rsync will fail even though all the
        ## preceeding tasks were sucessful
        if not os.path.exists(self.geno_concordance_out):
            os.system("touch %s"%self.geno_concordance_out)
        ## create the base directory if it does not exist
        if not os.path.exists(self.base_dir):
            os.makedirs(self.base_dir)

        if not os.path.exists(self.script_dir):
            os.makedirs(self.script_dir)
        ## refer to rsync man for more on options, its all the options in archive mode less preserving the permissions
        ## the o will not take effect unless run by a superuser and the g will default to your default group, i have removed
        ## the o option for safety. 
        if self.sample_type.upper() != 'GENOME': ## There is another additional capture specificity file 
            cmd = ("rsync --timeout=25000 -vrltgD "
                   "{script_dir} {log_dir} {recal_bam} {recal_bam_index} {annotated_vcf_gz} "
                   "{annotated_vcf_gz_index} {g_vcf_gz} {g_vcf_gz_index} "
                   "{cvg_binned} {gq_binned} {alignment_out} {cvg_out} "
                   "{cvg_ccds_out} {cvg_X_out} {cvg_Y_out} {cvg_cs_out} {dup_out} "
                   "{variant_call_out} {geno_concordance_out} {contamination_out} {fastq_loc} {base_dir}"
                   ).format(**self.__dict__)
        else:
            cmd = ("rsync --timeout=25000 -vrltgoD "
                   "{script_dir} {log_dir} {recal_bam} {recal_bam_index} {annotated_vcf_gz} "
                   "{annotated_vcf_gz_index} {g_vcf_gz} {g_vcf_gz_index} "
                   "{cvg_binned} {gq_binned} {alignment_out} {cvg_out} "
                   "{cvg_ccds_out} {cvg_X_out} {cvg_Y_out} {dup_out} "
                   "{variant_call_out} {geno_concordance_out} {contamination_out} {fastq_loc} {base_dir}"
                  ).format(**self.__dict__)
            if os.path.exists(self.original_vcf_gz):
                cmd = ("rsync --timeout=25000 -vrltgoD "
                   "{script_dir} {log_dir} {recal_bam} {recal_bam_index} {annotated_vcf_gz} "
                   "{annotated_vcf_gz_index} {g_vcf_gz} {g_vcf_gz_index} "
                   "{cvg_binned} {gq_binned} {alignment_out} {cvg_out} "
                   "{cvg_ccds_out} {cvg_X_out} {cvg_Y_out} {dup_out} "
                    "{variant_call_out} {geno_concordance_out} {contamination_out} {fastq_loc} {original_vcf_gz} {original_vcf_gz_index} {base_dir} "
                  ).format(**self.__dict__)
        try:
            ## Insert start time to database
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=self.pseudo_prepid,
                        pipeline_step_id=self.pipeline_step_id))
            db.commit()
            db.close()

            with open(self.script,'w') as o:
                o.write(cmd + "\n")
            DEVNULL = open(os.devnull, 'w')

            ## Run the rsync
            subprocess.check_call(shlex.split(cmd), stdout=DEVNULL,stderr=subprocess.STDOUT,close_fds=True)
            ## Change folder permissions of base directory
            subprocess.check_output("chgrp bioinfo {0}".format(self.base_dir),shell=True)
            subprocess.check_output("chmod -R 775 {0}".format(self.base_dir),shell=True)

            """Original dragen BAM could be technically deleted earlier after the
            realigned BAM has been created on scratch space but it is safer to
            delete after the final realigned, recalculated BAM has been archived
            since our scratch space has failed in the past."""
            rm_cmd = ['rm',self.bam,self.bam_index]
            rm_folder_cmd = ['rm','-rf',self.scratch_dir]
            #subprocess.call(rm_cmd)
            subprocess.call(rm_folder_cmd)

            DBID,prepID = getDBIDMaxPrepID(self.pseudo_prepid)
            userID = getUserID()
            ## Update the AlignSeqFile loc to the final archive location
            location = "{0}/{1}".format(config().base,self.sample_type.upper())
            update_alignseqfile(location,self.sample_name,self.pseudo_prepid)

            ## Insert to dragen status and update finish time in database
            db = get_connection("seqdb")
            cur = db.cursor()
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
            db.close()
            DEVNULL.close()
        finally:
            if db.open:
                db.close()

    def requires(self):
        yield self.clone(AnnotateVCF)
        if (self.sample_name.upper().startswith('SRR')):
            yield self.clone(SubsetVCF)
        yield self.clone(GQBinning)
        yield self.clone(CvgBinning)
        yield self.clone(UpdateSeqdbMetrics)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)


class SQLTarget(luigi.Target):
    """ A luigi target class describing verification of the entries in the database
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

############################################ POST GATK MODULES ##############################################
#########       
#########                              ADD COMMENTS HERE ........                                              
#########
######### 
#############################################################################################################

def get_pipeline_step_id(step_name,db_name):
    """
    Get the pipeline stepid from seqdb, will try 4 more
    times to execute the query at set intervals incase 
    the first attempt fails. Raises an approriate error
    if all the attempts failed.

    Args : 
         step_name : String ; name of the luigi Task Class name
         db_name : String ; Name of the database to query, either seqdb or
                   testdb

    Returns :
         The pipeline stepid.
    """

    tries = 5
    for i in range(tries):
        try:
            db = get_connection(db_name)
            cur = db.cursor()
            cur.execute(GET_PIPELINE_STEP_ID.format(
                step_name=step_name))
            pipeline_step_id = cur.fetchone()[0]
        except MySQLdb.Error, e:
            if i == tries - 1:
                raise Exception("ERROR %d IN CONNECTION: %s" % (e.args[0], e.args[1]))
            else:
                time.sleep(60) ## Wait a minute before trying again
                continue
        finally:
            if db.open:
                db.close()
            return pipeline_step_id ## Pass control back 

def run_shellcmd(db_name,pipeline_step_id,pseudo_prepid,cmd,sample_name,
                 task_name):
    """
    Runs a shell command, is flanked by database update calls
    to ensure logging the time and to check for task completion.
    If the database connection fails, will try 4 more attempts

    Args :
         db_name : String ; Name of the database to query, either seqdb
                   or testdb
         pipeline_step_id : String; Id for the luigi task in the database
         pseudo_prepid : Int; An indentifier for the sample in the database
         cmd : List; The command(s) to be run
         sample_name : String; CHGVID
         task_name : The luigi task name

    Returns :
         Does not return anything
    """

    tries = 5
    for i in range(tries):
        try:
            db = get_connection(db_name)
            cur = db.cursor()
            cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                        pseudo_prepid=pseudo_prepid,
                        pipeline_step_id=pipeline_step_id))
            db.commit()
            db.close()

            for cmd_j in cmd:
                proc = subprocess.Popen(cmd_j,shell=True)
                proc.wait()
                if proc.returncode: ## Non zero return code
                    raise subprocess.CalledProcessError(proc.returncode,cmd_j)

            db = get_connection(db_name)
            cur = db.cursor()
            update_dragen_status(pseudo_prepid,sample_name,task_name,cur)
            cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                        pseudo_prepid=pseudo_prepid,
                        pipeline_step_id=pipeline_step_id))
            db.commit()
            db.close()
        except MySQLdb.Error, e:
                if i == tries - 1:
                    raise Exception("ERROR %d IN CONNECTION: %s" % (e.args[0], e.args[1]))
                else:
                    time.sleep(60) ## Wait a minute before trying again
                    continue
        finally:
            if db.open:
                db.close()
            break ## Exit the loop 


def update_dragen_status(pseudo_prepid,sample_name,task_name,cur):
    """
    Updates dragen_statusT table with the pipeline step name 
    Note : This function calls getDBIDMaxPrepID and getUserID
    to get the DBID,prepID and userID 
    
    Args : 
         pseudo_prepid : Int; An indentifier for the sample in the database
         sample_name : String; CHGVID, a string identifier for the sequecing sample
         task_name : String; Luigi task name, will be used to update the status field
         cur : A MySQLdb cursor instance  
    Returns :
         Does not return anything
    """

    DBID,prepID = getDBIDMaxPrepID(pseudo_prepid)
    userID = getUserID()
    cur.execute(INSERT_DRAGEN_STATUS.format(
        sample_name=sample_name,status="Dragen "+task_name,DBID=DBID,
        prepID=prepID,pseudoPrepID=pseudo_prepid,userID=userID))

    
def get_productionvcf(pseudo_prepid,sample_name,sample_type):
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
                           sample_name,sample_type,pseudo_prepid)
                      )
    tries = 5
    for i in range(tries):
        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(query_statement)
            db_val = cur.fetchall()
            if len(db_val) == 0: ## If the query returned no results
                return None ## Pass control back, note the finally clause is still executed 
            if len(db_val) > 1:
                warnings.warn("More than 1 entry , warning :" 
                              "duplicate pseudo_prepids !")
            alignseqfileloc = db_val[-1][0] ## Get the last result
            vcf_loc = os.path.join(alignseqfileloc,'combined')
            temp_vcf = glob.glob(os.path.join(vcf_loc,
                                              '{0}.analysisReady.annotated.vcf.*'.format(sample_name)))
            if len(temp_vcf) == 0:
                return None 
            vcf = temp_vcf[0]
            if not os.path.isfile(vcf):
                vcf = os.path.join(vcf_loc,"{0}.analysisReady.annotated.vcf".format(sample_name))
                if not os.path.isfile(vcf):
                    return None
                else:
                    return vcf
            else:
                return vcf 
        except MySQLdb.Error, e:
            if i == tries - 1:
                raise Exception("ERROR %d IN CONNECTION: %s" % (e.args[0], e.args[1]))
            else:
                time.sleep(60) ## Wait a minute before trying again
                continue 
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

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    ## System Parameters
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self,*args,**kwargs):
        super(CreateGenomeBed,self).__init__(*args,**kwargs)
        self.seqtype_dir = os.path.join(self.scratch,self.sample_type.upper())
        self.sample_dir = os.path.join(self.seqtype_dir,self.sample_name+'.'+self.pseudo_prepid)
        kwargs["sample_name"] = self.sample_name
        super(CreateGenomeBed, self).__init__(*args, **kwargs)

        self.output_dir = os.path.join(self.sample_dir,'cvg_binned')
        self.genomecov_bed = os.path.join(self.sample_dir,self.sample_name+'.genomecvg.bed')
        self.recal_bam = os.path.join(self.sample_dir,self.sample_name+'.{0}.realn.recal.bam'.format(self.pseudo_prepid))
        self.log_dir = os.path.join(self.sample_dir,'logs')
        self.log_file = os.path.join(self.log_dir,self.sample_name+'.{0}.genomecov.log'.format(self.pseudo_prepid))
        self.genomecov_cmd = "{0} genomecov -bga -ibam {1} 1> {2} 2> {3}".format(
                              config().bedtools,self.recal_bam,
                              self.genomecov_bed,self.log_file)
        self.pipeline_step_id = get_pipeline_step_id(
            self.__class__.__name__,"seqdb")

    def requires(self):
        """
        The dependency is the completion of the HaplotypeCaller task
        Note : Can use PrintReads here to speed things up, but may
        face IO bottle necks --> Future test
        """

        return self.clone(HaplotypeCaller)

    def output(self):
        """
        Output is the genomecov bed file
        """
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
            pipeline_step_id=self.pipeline_step_id)

    def work(self):
        """
        Run the bedtools cmd
        """
        run_shellcmd("seqdb",self.pipeline_step_id,self.pseudo_prepid,
                     [self.genomecov_cmd],self.sample_name,self.__class__.__name__)

class CalculateTargetCoverage(SGEJobTask):
    """ Run GATK's calculate target coverage for getting variables for CNV 
    analysis
    """
    
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    ## System Parameters
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"


    def __init__(self,*args,**kwargs):
        super(CalculateTargetCoverage,self).__init__(*args,**kwargs)
        self.seqtype_dir = "{0}/{1}".format(self.scratch,self.sample_type.upper())
        self.sample_dir = os.path.join(self.seqtype_dir,self.sample_name+'.'+self.pseudo_prepid)
        self.padded_target_cvg_file = "{0}/{1}.{2}/{1}.{2}.padded.target.cov.tsv".format(self.seqtype_dir,self.sample_name,self.pseudo_prepid)
        self.recal_bam = "{0}/{1}.{2}/{1}.{2}.realn.recal.bam".format(self.seqtype_dir,self.sample_name,self.pseudo_prepid)
        self.script = "{0}/{1}.{2}/scripts/{1}.{2}.{3}.sh".format(
            self.seqtype_dir,self.sample_name,self.pseudo_prepid,self.__class__.__name__)
        kwargs["sample_name"] = self.sample_name
        super(CalculateTargetCoverage, self).__init__(*args, **kwargs)

        self.calc_target_cov_cmd = ("{0} -jar {1} CalculateTargetCoverage "
            "--output {2} "
            "--input {3} "
            "--reference {4} "
            "--cohortName {5}.{6} "
            "--interval_padding {7} "
            "-L {8} "
            "--keepduplicatereads true "
            "--disable_all_read_filters false "
            "--interval_set_rule UNION "
            "--targetInformationColumns FULL "
            "--readValidationStringency SILENT "
            ).format(config().java,
                    config().gatk4,
                    self.padded_target_cvg_file,
                    self.recal_bam,
                    config().ref,
                    self.sample_name,
                    self.pseudo_prepid,
                    config().cnv_target_padding,
                    self.capture_kit_bed)

        self.pipeline_step_id = get_pipeline_step_id(
            self.__class__.__name__,"seqdb")

    def requires(self):
        if self.sample_type.lower() == 'exome':
            return self.clone(HaplotypeCaller)
        else:
            raise Exception, "Sample type: %s not supported in this module" % self.sample_type
        """
        elif self.sample_type.lower() == 'genome':
            return self.clone(HaplotypeCaller)
        """
    def output(self):
        self.pipeline_step_id = get_pipeline_step_id(
            self.__class__.__name__,"seqdb")
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

    def work(self):

        with open(self.script,'w') as o:
            o.write(self.calc_target_cov_cmd + "\n")

        """
        Run the GATK's CalculateTargetCoverage cmd for CNV calling
        """
        run_shellcmd("seqdb",self.pipeline_step_id,self.pseudo_prepid,
                     [self.calc_target_cov_cmd],self.sample_name,self.__class__.__name__)

class NormalizeSomaticReadCounts(SGEJobTask):
    """ Normalization step in the CNV pipeline
    """
    
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    ## System Parameters
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self,*args,**kwargs):
        super(NormalizeSomaticReadCounts,self).__init__(*args,**kwargs)
        self.seqtype_dir = "{0}/{1}".format(self.scratch,self.sample_type.upper())
        self.padded_target_cvg_file = "{0}/{1}.{2}/{1}.{2}.padded.target.cov.tsv".format(self.seqtype_dir,self.sample_name,self.pseudo_prepid)
        self.preTangentNormalized = "{0}/{1}.{2}/{1}.{2}.ptn.tsv".format(self.seqtype_dir,self.sample_name,self.pseudo_prepid)
        self.tangentNormalized = "{0}/{1}.{2}/{1}.{2}.tn.tsv".format(self.seqtype_dir,self.sample_name,self.pseudo_prepid)
        kwargs["sample_name"] = self.sample_name
        super(NormalizeSomaticReadCounts, self).__init__(*args, **kwargs)

        capture_kit = getCaptureKit(self.pseudo_prepid)
        panelOfNormals = getPanelOfNormals(capture_kit,self.sample_type)
        panelOfNormals = '/nfs/seqscratch11/normal100bp.pon'
        self.normalize_cmd = ("{java} -jar {gatk4} NormalizeSomaticReadCounts "
            "-I {padded_target_cvg_file} "
            "-PON {panelOfNormals} "
            "-PTN {preTangentNormalized} "
            "-TN {tangentNormalized}"
            ).format(java=config().java,gatk4=config().gatk4,
                    sample_name=self.sample_name,
                    pseudo_prepid=self.pseudo_prepid,
                    preTangentNormalized=self.preTangentNormalized,
                    tangentNormalized=self.tangentNormalized,
                    padded_target_cvg_file=self.padded_target_cvg_file,
                    panelOfNormals=panelOfNormals)

        self.pipeline_step_id = get_pipeline_step_id(
            self.__class__.__name__,"seqdb")

    def requires(self):
        return self.clone(CalculateTargetCoverage)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

    def work(self):
        print self.normalize_cmd
        """
        Run the GATK's NormalizeSomaticReadCounts cmd for CNV calling
        """
        run_shellcmd("seqdb",self.pipeline_step_id,self.pseudo_prepid,
                     [self.normalize_cmd],self.sample_name,self.__class__.__name__)

class PerformSegmentation(SGEJobTask):
    """ Perform segmentation of normalized copy number ratios
    """
    
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    ## System Parameters
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self,*args,**kwargs):
        super(PerformSegmentation,self).__init__(*args,**kwargs)
        self.seqtype_dir = "{0}/{1}".format(self.scratch,self.sample_type.upper())
        self.padded_target_cvg_file = "{0}/{1}.{2}/{1}.{2}.padded.target.cov.tsv".format(self.seqtype_dir,self.sample_name,self.pseudo_prepid)
        self.preTangentNormalized = "{0}/{1}.{2}/{1}.{2}.ptn.tsv".format(self.seqtype_dir,self.sample_name,self.pseudo_prepid)
        self.tangentNormalized = "{0}/{1}.{2}/{1}.{2}.tn.tsv".format(self.seqtype_dir,self.sample_name,self.pseudo_prepid)
        self.segments = "{0}/{1}.{2}/{1}.{2}.seg".format(self.seqtype_dir,self.sample_name,self.pseudo_prepid)
        kwargs["sample_name"] = self.sample_name
        super(PerformSegmentation, self).__init__(*args, **kwargs)

        capture_kit = getCaptureKit(self.pseudo_prepid)
        panelOfNormals = getPanelOfNormals(capture_kit,self.sample_type)
        panelOfNormals = "/nfs/seqscratch11/normal100bp.pon"
        self.seg_cmd = ("{java} -jar {gatk4} PerformSegmentation "
            "-TN {tangentNormalized} "
            "-O {segments}"
            ).format(java=config().java,gatk4=config().gatk4,
                    tangentNormalized=self.tangentNormalized,
                    segments=self.segments,
                    panelOfNormals=panelOfNormals)

        self.pipeline_step_id = get_pipeline_step_id(
            self.__class__.__name__,"seqdb")

    def requires(self):
        return self.clone(NormalizeSomaticReadCounts)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

    def work(self):
        print self.seg_cmd
        run_shellcmd("seqdb",self.pipeline_step_id,self.pseudo_prepid,
                     [self.seg_cmd],self.sample_name,self.__class__.__name__)

class CvgBinning(SGEJobTask):
    """
    Task to run the binning script for Coverage
    """

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self,*args,**kwargs):

        super(CvgBinning,self).__init__(*args,**kwargs)
        self.seqtype_dir = os.path.join(self.scratch,self.sample_type.upper())
        self.sample_dir = os.path.join(self.seqtype_dir,self.sample_name+'.'+self.pseudo_prepid)
        kwargs["sample_name"] = self.sample_name
        super(CvgBinning,self).__init__(*args, **kwargs)

        self.output_dir = os.path.join(self.sample_dir,'cvg_binned')
        self.genomecov_bed = os.path.join(self.sample_dir,self.sample_name+'.genomecvg.bed')
        self.log_dir = os.path.join(self.sample_dir,'logs')
        self.log_file = os.path.join(self.log_dir,self.sample_name+'.{0}.cvgbinning.log'.format(self.pseudo_prepid))
        self.script = os.path.join(self.sample_dir,"scripts/{0}.{1}.cvg.binning.sh".format(self.sample_name,self.pseudo_prepid))
        try: ## This syntax is needed to avoid race conditions in certain cases  
            if not os.path.isdir(self.output_dir): ## Recursively create the directory if it doesnt exist
                os.makedirs(self.output_dir)
        except OSError,e:
            if e.errno != 17:
                raise Exception("Problem Creating Directory : {0}".format(e))
            pass


        self.human_chromosomes = []
        self.human_chromosomes.extend(range(1, 23))
        self.human_chromosomes = [str(x) for x in self.human_chromosomes]
        self.human_chromosomes.extend(['X', 'Y','MT'])

        self.binning_cmd = "{0} {1} 1000 {2} {3} {4} >& {5}".format(config().pypy,config().coverage_binner,
                                                                    self.sample_name+'.'+self.pseudo_prepid,self.genomecov_bed,
                                                                    self.output_dir,self.log_file)
        self.rm_cmd = "rm {0}".format(self.genomecov_bed)
        self.pipeline_step_id = get_pipeline_step_id(
            self.__class__.__name__,"seqdb")

    def requires(self):
        """
        Dependency is the completion of the CreateGenomeBed Task
        """

        return self.clone(CreateGenomeBed)

    def output(self):
        """
        Output from this task are 23 files for each chromosome
        Also check for deletion of the genomecov file
        """

        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

    def work(self):
        """ Run the binning script
        """
        with open(self.script,'w') as OUT:
            OUT.write(self.binning_cmd+'\n')
        run_shellcmd("seqdb",self.pipeline_step_id,self.pseudo_prepid,
                     [self.binning_cmd,self.rm_cmd],self.sample_name,
                     self.__class__.__name__)

class GQBinning(SGEJobTask):
    """
    Task to run the binning script for GQ
    """
    
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"  
    
    def __init__(self,*args,**kwargs):
        super(GQBinning,self).__init__(*args,**kwargs)
        self.seqtype_dir = os.path.join(self.scratch,self.sample_type.upper())
        self.sample_dir = os.path.join(self.seqtype_dir,self.sample_name+'.'+self.pseudo_prepid)
        kwargs["sample_name"] = self.sample_name
        super(GQBinning,self).__init__(*args, **kwargs)
        
        self.output_dir = os.path.join(self.sample_dir,'gq_binned')
        self.gvcf = os.path.join(self.sample_dir,self.sample_name+'.{0}.g.vcf.gz'.format(self.pseudo_prepid))
        self.log_dir = os.path.join(self.sample_dir,'logs')
        self.log_file = os.path.join(self.log_dir,self.sample_name+'.{0}.gqbinning.log'.format(self.pseudo_prepid))
        self.script = os.path.join(self.sample_dir,"scripts/{0}.{1}.gqbinning.sh".format(self.sample_name,self.pseudo_prepid))
        try:
            if not os.path.isdir(self.output_dir): ## Recursively create the directory if it doesnt exist
                os.makedirs(self.output_dir)
        except OSError,e:
            if e.errno != 17:
                raise Exception("Problem Creating Directory : {0}".format(e))
            pass
        
    
        self.human_chromosomes = []
        self.human_chromosomes.extend(range(1, 23))
        self.human_chromosomes = [str(x) for x in self.human_chromosomes]
        self.human_chromosomes.extend(['X', 'Y','MT'])
                       

        self.binning_cmd = """ {0} {1} 10000 {2} {3} {4} >& {5} """.format(
            config().pypy,config().gq_binner,self.gvcf,self.sample_name+'.'+self.pseudo_prepid,self.output_dir,self.log_file)
        
        self.pipeline_step_id = get_pipeline_step_id(
            self.__class__.__name__,"seqdb")
        
    def requires(self):
        """
        Dependency is the completion of the HaplotypeCaller Task
        """
        
        return self.clone(HaplotypeCaller)
        
    def output(self):
        """
        Output from this task are 23 files, one for each chromosome
        """
        
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

    def work(self):
        with open(self.script,'w') as OUT:
            OUT.write(self.binning_cmd+'\n')
        run_shellcmd("seqdb",self.pipeline_step_id,self.pseudo_prepid,
                     [self.binning_cmd],self.sample_name,self.__class__.__name__)

class AlignmentMetrics(SGEJobTask):
    """ Run Picard AlignmentSummaryMetrics """

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self,*args,**kwargs):
        super(AlignmentMetrics,self).__init__(*args,**kwargs)
        self.seqtype_dir = os.path.join(self.scratch,self.sample_type.upper())
        self.sample_dir = os.path.join(self.seqtype_dir,self.sample_name+'.'+self.pseudo_prepid)
        kwargs["sample_name"] = self.sample_name
        super(AlignmentMetrics, self).__init__(*args, **kwargs)

        self.recal_bam = os.path.join(self.sample_dir,self.sample_name+'.{0}.realn.recal.bam'.format(self.pseudo_prepid))
        self.output_file_raw = os.path.join(self.sample_dir,self.sample_name+'.alignment.metrics.raw.txt')
        self.output_file = os.path.join(self.sample_dir,self.sample_name+'.{0}.alignment.metrics.txt'.format(self.pseudo_prepid))
        self.script = os.path.join(self.sample_dir,"scripts/{0}.{1}.alignment.metrics.sh".format(self.sample_name,self.pseudo_prepid))
        self.log_dir = os.path.join(self.sample_dir,'logs')
        self.log_file = os.path.join(self.log_dir,self.sample_name+'.{0}.alignment.metrics.log'.format(self.pseudo_prepid))
        self.cmd = "{0} -XX:ParallelGCThreads=4 -jar {1} CollectAlignmentSummaryMetrics TMP_DIR={2} VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE={3} INPUT={4} OUTPUT={5} >& {6}".format(config().java,config().picard,self.sample_dir,config().ref,self.recal_bam,self.output_file_raw,self.log_file)
        self.parser_cmd = """cat {0} | grep -v "^#" | awk -f {1} > {2}""".format(self.output_file_raw,config().transpose_awk,self.output_file)
        self.pipeline_step_id = get_pipeline_step_id(
            self.__class__.__name__,"seqdb")

    def output(self):
        """Check whether the output file is present"""
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

    def requires(self):
        """The dependencies for this task is simply the existence of the bam file
        from dragen with duplicates removed"""
        return self.clone(HaplotypeCaller)

    def work(self):
        self.remove_cmd = "rm {0}".format(self.output_file_raw)
        with open(self.script,'w') as OUT:
            OUT.write(self.cmd+'\n')
            OUT.write(self.parser_cmd+'\n')
            OUT.write(self.remove_cmd+'\n')
        run_shellcmd("seqdb",self.pipeline_step_id,self.pseudo_prepid,
                     [self.cmd,self.parser_cmd,self.remove_cmd],self.sample_name,
                     self.__class__.__name__)

class RunCvgMetrics(SGEJobTask):
    """ Task to get the coverage metrics"""

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self,*args,**kwargs):
        super(RunCvgMetrics,self).__init__(*args,**kwargs)
        ## Define on the fly parameters
        self.seqtype_dir = os.path.join(self.scratch,self.sample_type.upper())
        self.sample_dir = os.path.join(self.seqtype_dir,self.sample_name+'.'+self.pseudo_prepid)
        kwargs["sample_name"] = self.sample_name
        super(RunCvgMetrics, self).__init__(*args, **kwargs)

        self.recal_bam = os.path.join(self.sample_dir,self.sample_name+'.{0}.realn.recal.bam'.format(self.pseudo_prepid))
        self.output_file = os.path.join(self.sample_dir,self.sample_name + ".{0}.cvg.metrics.txt".format(self.pseudo_prepid))
        self.raw_output_file = os.path.join(self.sample_dir,self.sample_name + ".cvg.metrics.raw.txt")
        self.output_file_ccds = os.path.join(self.sample_dir,self.sample_name + ".{0}.cvg.metrics.ccds.txt".format(self.pseudo_prepid))
        self.raw_output_file_ccds = os.path.join(self.sample_dir,self.sample_name + ".cvg.metrics.ccds.raw.txt")
        self.output_file_X = os.path.join(self.sample_dir,self.sample_name+ ".{0}.cvg.metrics.X.txt".format(self.pseudo_prepid))
        self.raw_output_file_X = os.path.join(self.sample_dir,self.sample_name + ".cvg.metrics.X.raw.txt")
        self.output_file_Y = os.path.join(self.sample_dir,self.sample_name+ ".{0}.cvg.metrics.Y.txt".format(self.pseudo_prepid))
        self.raw_output_file_Y = os.path.join(self.sample_dir,self.sample_name + ".cvg.metrics.Y.raw.txt")
        self.raw_output_file_cs = os.path.join(self.sample_dir,self.sample_name + ".cvg.metrics.cs.raw.txt")
        self.output_file_cs = os.path.join(self.sample_dir,self.sample_name + ".{0}.cvg.metrics.cs.txt".format(self.pseudo_prepid))
        self.target_file = config().target_file ## Targets here refers to ccds 
        self.log_dir = os.path.join(self.sample_dir,'logs')
        self.script = os.path.join(self.sample_dir,"scripts/{0}.{1}.cvg.metrics.sh".format(self.sample_name,self.pseudo_prepid))
        self.log_file = os.path.join(self.log_dir,self.sample_name+'.{0}.cvg.metrics.log'.format(self.pseudo_prepid))

        self.DBID,self.prepID = getDBIDMaxPrepID(self.pseudo_prepid)


        if self.sample_type.upper() == "GENOME":
            ## Define shell commands to be run
            self.cvg_cmd = """{0} -Xmx{1}g -XX:ParallelGCThreads=4 -jar {2} CollectWgsMetrics VALIDATION_STRINGENCY=LENIENT R={3} I={4} INTERVALS={5} O={6} MQ=20 Q=10 >>{6}"""
            ## Run on the ccds regions only
            self.cvg_cmd1 = self.cvg_cmd.format(config().java,config().max_mem,config().picard,config().ref,self.recal_bam,config().target_file,self.raw_output_file_ccds,self.log_file)
            ## Run across the genome
            self.cvg_cmd2 = """{0} -Xmx{1}g -XX:ParallelGCThreads=4 -jar {2} CollectWgsMetrics R={3} I={4} O={5} MQ=20 Q=10 >> {6}""".format(config().java,config().max_mem,config().picard,config().ref,self.recal_bam,self.raw_output_file,self.log_file)
            ## Run on X and Y Chromosomes only (across all regions there not just ccds)
            self.cvg_cmd3 = self.cvg_cmd.format(config().java,config().max_mem,config().picard,config().ref,self.recal_bam,config().target_file_X,self.raw_output_file_X,self.log_file)
            self.cvg_cmd4 = self.cvg_cmd.format(config().java,config().max_mem,config().picard,config().ref,self.recal_bam,config().target_file_Y,self.raw_output_file_Y,self.log_file)
        else:
            db = get_connection("seqdb")
            cur = db.cursor()
            query = """SELECT region_file_lsrc FROM captureKit INNER JOIN prepT ON prepT_name=prepT.exomeKit WHERE prepID = {0} AND (chr = 'X' OR chr = 'Y')""".format(self.prepID)
            cur.execute(query)
            db_val = cur.fetchall()
            self.capture_file_X = db_val[0][0]
            self.capture_file_Y = db_val[1][0]
            db.close()

            self.calc_capture_specificity = 1 ## We need to run an additional picard module for calculating capture specificity
            self.cvg_cmd = "{0} -Xmx{1}g -XX:ParallelGCThreads=4 -jar {2} -T DepthOfCoverage -mbq 10 -mmq 20 --omitIntervalStatistics --omitDepthOutputAtEachBase --omitLocusTable -R {3} -I {4} -ct 5 -ct 10 -ct 15 -ct 20 -L {5} -o {6} >> {7}"
            ## Run across ccds regions
            self.cvg_cmd1 = self.cvg_cmd.format(config().java,config().max_mem,config().gatk,config().ref,self.recal_bam,config().target_file,self.raw_output_file_ccds,self.log_file)
            ## Run across capture kit regions
            self.cvg_cmd2 = self.cvg_cmd.format(config().java,config().max_mem,config().gatk,config().ref,self.recal_bam,self.capture_kit_bed,self.raw_output_file,self.log_file)
            ## Run on X and Y Chromosomes only (across all regions there not just ccds)
            self.cvg_cmd3 = self.cvg_cmd.format(config().java,config().max_mem,config().gatk,config().ref,self.recal_bam,self.capture_file_X,self.raw_output_file_X,self.log_file)
            self.cvg_cmd4 = self.cvg_cmd.format(config().java,config().max_mem,config().gatk,config().ref,self.recal_bam,self.capture_file_Y,self.raw_output_file_Y,self.log_file)
            ## Run PicardHsMetrics for Capture Specificity
            self.output_bait_file = os.path.join(self.sample_dir,'targets.interval')
            self.convert_cmd = "{0} -jar {1} BedToIntervalList I={2} SD={3} OUTPUT={4} >> {5}".format(config().java,config().picard,self.capture_kit_bed,config().seqdict_file,self.output_bait_file,self.log_file)
            self.cvg_cmd5 = """{0} -Xmx{1}g -XX:ParallelGCThreads=4 -jar {2} CollectHsMetrics BI={3} TI={4} VALIDATION_STRINGENCY=SILENT METRIC_ACCUMULATION_LEVEL=ALL_READS I={5} O={6} MQ=20 Q=10 >> {7}""".format(config().java,config().max_mem,config().picard,self.output_bait_file,self.output_bait_file,self.recal_bam,self.raw_output_file_cs,self.log_file)
            ## DepthOfCoverage output has an extra suffix at the end
            self.raw_output_file_ccds = self.raw_output_file_ccds+'.sample_summary'
            self.raw_output_file = self.raw_output_file+'.sample_summary'
            self.raw_output_file_X = self.raw_output_file_X+'.sample_summary'
            self.raw_output_file_Y = self.raw_output_file_Y+'.sample_summary'

        ## An intermediate parsed file is created for extracting info for db update later, this remains the same for exomes and genomes
        self.parser_cmd = """cat {0} | grep -v "^#" | awk -f {1} > {2}"""
        self.parser_cmd1 = self.parser_cmd.format(self.raw_output_file_ccds,config().transpose_awk,self.output_file_ccds)
        self.parser_cmd2 = self.parser_cmd.format(self.raw_output_file,config().transpose_awk,self.output_file)
        self.parser_cmd3 = self.parser_cmd.format(self.raw_output_file_X,config().transpose_awk,self.output_file_X)
        self.parser_cmd4 = self.parser_cmd.format(self.raw_output_file_Y,config().transpose_awk,self.output_file_Y)
        if self.sample_type.upper() != "GENOME": ## i.e. Exome or CustomCapture
            self.parser_cmd5 = self.parser_cmd.format(self.raw_output_file_cs,config().transpose_awk,self.output_file_cs)
        self.pipeline_step_id = get_pipeline_step_id(
            self.__class__.__name__,"seqdb")

    def output(self):
        """The output produced by this task"""
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

    def requires(self):
        """The dependency for this task is the PrintReads task"""
        yield self.clone(HaplotypeCaller)

    def work(self):
        """Run Picard CalculateHsMetrics,DepthOfCoverage or CollectWgsMetrics
        An optimization for the future is to split this up into multiple tasks since they are all independent
        of each other will speed things up significantly, IO can be an issue"""

        with open(self.script,'w') as OUT:
            OUT.write(self.cvg_cmd1+'\n')
            OUT.write(self.parser_cmd1+'\n')
            OUT.write(self.cvg_cmd2+'\n')
            OUT.write(self.parser_cmd2+'\n')
            OUT.write(self.cvg_cmd3+'\n')
            OUT.write(self.parser_cmd3+'\n')
            OUT.write(self.cvg_cmd4+'\n')
            OUT.write(self.parser_cmd4+'\n')

        commands_run = [self.cvg_cmd1,self.parser_cmd1,self.cvg_cmd2,
                      self.parser_cmd2,self.cvg_cmd3,self.parser_cmd3,
                      self.cvg_cmd4,self.parser_cmd4]

        if self.sample_type.upper()!='GENOME': ## Append additional commands for capture specificity
            with open(self.script,'w') as OUT:
                OUT.write(self.convert_cmd+'\n')
                OUT.write(self.cvg_cmd5+'\n')
                OUT.write(self.parser_cmd5+'\n')
            commands_run.extend([self.convert_cmd,self.cvg_cmd5,self.parser_cmd5])

        ## Add the remove commands for clearing temp files
        #remove_cmd = "rm {0}/*.cvg.metrics*raw*".format(self.sample_dir)
        #commands_run.append(remove_cmd)
        run_shellcmd("seqdb",self.pipeline_step_id,self.pseudo_prepid,commands_run,
                     self.sample_name,self.__class__.__name__)

class DuplicateMetrics(SGEJobTask):
    """ Parse Duplicate Metrics from dragen logs """

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self,*args,**kwargs):
        super(DuplicateMetrics,self).__init__(*args,**kwargs)
        self.seqtype_dir_seqscratch = os.path.join(self.scratch,self.sample_type.upper())
        self.sample_dir_seqscratch = os.path.join(self.seqtype_dir_seqscratch,self.sample_name+'.'+self.pseudo_prepid)
        self.seqtype_dir = os.path.join(self.scratch,self.sample_type.upper())
        self.sample_dir = os.path.join(self.seqtype_dir,self.sample_name+'.'+self.pseudo_prepid)
        kwargs["sample_name"] = self.sample_name
        super(DuplicateMetrics, self).__init__(*args, **kwargs)
        self.log_dir = os.path.join(self.sample_dir,'logs')
        self.log_dir_seqscratch = os.path.join(self.sample_dir_seqscratch,'logs')
        self.dragen_log = os.path.join(self.log_dir,self.sample_name+'.'+self.pseudo_prepid+'.dragen.out')
        self.dragen_log_seqscratch = os.path.join(self.log_dir_seqscratch,self.sample_name+'.'+self.pseudo_prepid+'.dragen.out')
        self.old_dragen_log = os.path.join(self.log_dir,self.sample_name+'.'+self.pseudo_prepid+'.dragen.log')
        self.old_dragen_log_seqscratch = os.path.join(self.log_dir_seqscratch,self.sample_name+'.'+self.pseudo_prepid+'.dragen.log')
        self.output_file = os.path.join(self.sample_dir_seqscratch,self.sample_name+'.{0}.duplicates.txt'.format(self.pseudo_prepid))
        self.super_old_dragen_log = os.path.join(self.log_dir,self.sample_name+'.out')
        self.log_file = os.path.join(self.log_dir_seqscratch,self.sample_name+'.{0}.dups.log'.format(self.pseudo_prepid))

        ## Check for the dragen log_file in either the seqscratch or the fastq16 directory
        ## Just use a wild card search ;-/ 
        temp = glob.glob(os.path.join(self.log_dir_seqscratch,'*.dragen.out'))
        if len(temp) > 0:
            self.dlog = temp[0]
        else:
            temp = glob.glob(os.path.join(self.log_dir_seqscratch,'*.dragen.log'))
            if len(temp) > 0:
                self.dlog = temp[0]
            else:
                temp = glob.glob(os.path.join(self.log_dir,'*.dragen.out'))
                if len(temp) > 0:
                    self.dlog = temp[0]
                else:
                    temp = glob.glob(os.path.join(self.log_dir,'*.dragen.log'))
                    if len(temp) > 0:
                        self.dlog = temp[0]
                    else:
                        temp = glob.glob(os.path.join(self.log_dir,'*.out'))
                        if len(temp) > 0:
                            self.dlog = temp[0]
                        else:
                            temp = glob.glob(os.path.join(self.log_dir,"{0}.log".format(self.sample_name)))
                            if len(temp) > 0:
                                self.dlog = temp[0]
                            else:
                                self.dlog = self.dragen_log

        if not os.path.exists(self.dlog):
            raise Exception("The dragen log file could not be found !")

        self.cmd = """grep 'duplicates marked' %s"""%(self.dlog)
        self.cmd2 = """grep 'Number of duplicate reads' %s|head -1"""%(self.dlog)

        self.pipeline_step_id = get_pipeline_step_id(
            self.__class__.__name__,"seqdb")

    def requires(self):
        """ Dependencies for this task is the existence of the log file """

        #return MyExtTask(self.dlog)
        return self.clone(CopyBam)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

    def work(self):
        ## The regular expression to catch the percentage duplicates in the grep string
        catch_dup = re.compile('.*\((.*)%\).*')
        catch_dup2 = re.compile('.*\[(.*)\].*')
        flag = 0
        ## Try to find a more elegant way to do this
        tries = 5
        for i in range(tries):
            try:
                db = get_connection("seqdb")
                cur = db.cursor()
                cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                    pseudo_prepid=self.pseudo_prepid,
                    pipeline_step_id=self.pipeline_step_id))
                db.commit()
                db.close()

                proc=subprocess.Popen(self.cmd,shell=True,stdout=subprocess.PIPE)
                proc.wait()
                if proc.returncode: ## Non zero return code
                    flag = 1 
                    #raise subprocess.CalledProcessError(proc.returncode,self.cmd)

                if flag == 0:
                    dragen_output = proc.stdout.read()
                    match = re.match(catch_dup,dragen_output)
                    if match:
                        perc_duplicates = match.group(1)
                        with open(self.output_file,'w') as OUT_HANDLE:
                            print >> OUT_HANDLE,perc_duplicates
                    else:
                        flag = 1

                if flag == 1: ## Either the match was unsucessful or the grep was unsucessful

                    proc=subprocess.Popen(self.cmd2,shell=True,stdout=subprocess.PIPE)
                    proc.wait()

                    if proc.returncode: ## Non zero return code
                        with open(self.log_file,'w') as LOG:
                            print >> LOG,"Error in : {0} {1}".format(self.cmd2,str(proc.returncode))
                            raise subprocess.CalledProcessError(proc.returncode,self.cmd)
                    dragen_output = proc.stdout.read().strip('\n')
                    match = re.match(catch_dup2,dragen_output)                   

                    if match:
                        with open(self.log_file,'w') as LOG:
                            print >> LOG,match.group(1)

                        prec_duplicates = match.group(1)
                        with open(self.output_file,'w') as OUT_HANDLE:
                            print >> OUT_HANDLE,match.group(1)
                    else:
                        with open(self.log_file,'w') as LOG:
                            raise Exception("Could not find duplicate metrics in dragen log!")
                            print >> LOG,"Could not find duplicate metrics in dragen log!"

                db = get_connection("seqdb")
                cur = db.cursor()
                update_dragen_status(self.pseudo_prepid,self.sample_name,
                                     self.__class__.__name__,cur)
                cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                    pseudo_prepid=self.pseudo_prepid,
                    pipeline_step_id=self.pipeline_step_id))
                db.commit()
                db.close()

            except MySQLdb.Error, e:
                if i == tries - 1:
                    raise Exception("ERROR %d IN CONNECTION: %s" % (e.args[0], e.args[1]))
                else:
                    time.sleep(60) ## Wait a minute before trying again
                    continue
            finally:
                if db.open:
                    db.close()
                break ## Exit the loop since we have completed the task 

class VariantCallingMetrics(SGEJobTask):
    """ Run the picard tool for qc evaluation of the variant calls"""
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self,*args,**kwargs):
        super(VariantCallingMetrics,self).__init__(*args,**kwargs)
        self.seqtype_dir = os.path.join(self.scratch,self.sample_type.upper())
        self.sample_dir = os.path.join(self.seqtype_dir,self.sample_name+'.'+self.pseudo_prepid)
        kwargs["sample_name"] = self.sample_name
        super(VariantCallingMetrics, self).__init__(*args, **kwargs)
        self.output_file_raw = os.path.join(self.sample_dir,self.sample_name+'.raw')
        self.output_file_raw1 = os.path.join(self.sample_dir,self.sample_name+".raw.variant_calling_summary_metrics")
        self.output_file_raw2 = os.path.join(self.sample_dir,self.sample_name+".raw.variant_calling_detail_metrics")
        self.annotated_vcf_gz = os.path.join(self.sample_dir,"{0}.{1}.analysisReady.annotated.vcf.gz".format(self.sample_name,self.pseudo_prepid))
        self.output_file = os.path.join(self.sample_dir,"{0}.{1}.variant_calling_summary_metrics".format(self.sample_name,self.pseudo_prepid))
        self.script = os.path.join(self.sample_dir,"scripts/{0}.{1}.variantcalling.metrics.sh".format(self.sample_name,self.pseudo_prepid))
        self.log_dir = os.path.join(self.sample_dir,'logs')
        self.log_file = os.path.join(self.log_dir,self.sample_name+'.{0}.variantcalling.metrics.log'.format(self.pseudo_prepid))
        self.cmd = "{0} -XX:ParallelGCThreads=4 -jar {1} CollectVariantCallingMetrics INPUT={2} OUTPUT={3} DBSNP={4} >& {5}".format(config().java,config().picard,self.annotated_vcf_gz,self.output_file_raw,config().dbSNP,self.log_file)
        self.parser_cmd = """cat {0} | grep -v "^#" | awk -f {1}  > {2}""".format(self.output_file_raw1,config().transpose_awk,self.output_file)
        self.pipeline_step_id = get_pipeline_step_id(
            self.__class__.__name__,"seqdb")

    def requires(self):
        """ The requirement for this task is the presence of the analysis ready
        vcf file from the GATK pipeline """

        return self.clone(AnnotateVCF)

    def output(self):
        """ The result from this task is the creation of the metrics file """

        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

    def work(self):
        self.remove_cmd1 = "rm {0}".format(self.output_file_raw1)
        self.remove_cmd2 = "rm {0}".format(self.output_file_raw2)
        with open(self.script,'w') as OUT:
            OUT.write(self.cmd+'\n')
            OUT.write(self.parser_cmd+'\n')
            OUT.write(self.remove_cmd1+'\n')
            OUT.write(self.remove_cmd2+'\n')

        run_shellcmd("seqdb",self.pipeline_step_id,self.pseudo_prepid,
                     [self.cmd,self.parser_cmd,self.remove_cmd1,
                      self.remove_cmd2],self.sample_name,self.__class__.__name__)

class GenotypeConcordance(SGEJobTask):
    """ Run the gatk tool for evaluation of the variant calls
    with older production pipeline variant calls """

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self,*args,**kwargs):
        super(GenotypeConcordance,self).__init__(*args,**kwargs)
        self.seqtype_dir = os.path.join(self.scratch,self.sample_type.upper())
        self.sample_dir = os.path.join(self.seqtype_dir,self.sample_name+'.'+self.pseudo_prepid)
        kwargs["sample_name"] = self.sample_name
        super(GenotypeConcordance, self).__init__(*args, **kwargs)
        self.output_file = os.path.join(self.sample_dir,
                                        "{0}.{1}.genotype_concordance_metrics"
                                        .format(self.sample_name,self.pseudo_prepid))
        self.annotated_vcf_gz = os.path.join(self.sample_dir,"{0}.{1}.analysisReady.annotated.vcf.gz".format(self.sample_name,self.pseudo_prepid))
        self.script = os.path.join(self.sample_dir,"scripts/{0}.{1}.genotype.concordance.metrics.sh".format(self.sample_name,self.pseudo_prepid))
        self.log_dir = os.path.join(self.sample_dir,'logs')
        self.log_file = os.path.join(self.log_dir,
                                     self.sample_name+
                                     '.{0}.genotypeconcoradance.metrics.log'.format(self.pseudo_prepid))
        self.concordance_cmd = "{0} -XX:ParallelGCThreads=4 -jar {1} -T GenotypeConcordance -R {2} -eval {3} -comp {4} -o {5} -U ALLOW_SEQ_DICT_INCOMPATIBILITY >& {6}"
        self.pipeline_step_id = get_pipeline_step_id(
            self.__class__.__name__,"seqdb")

    def requires(self):
        """ The requirement for this task is the presence of the analysis ready
        vcf file from the GATK pipeline """
        return self.clone(AnnotateVCF)

    def output(self):
        """ The result from this task is the creation of the metrics file """

        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

    def work(self):
        ## Try to find a more elegant way to do this later
        tries = 5
        for i in range(tries):
            try:
                db = get_connection("seqdb")
                cur = db.cursor()
                cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                    pseudo_prepid=self.pseudo_prepid,
                    pipeline_step_id=self.pipeline_step_id))
                db.commit()
                db.close()
                ## The actual task commands to be run
                production_vcf = get_productionvcf(self.pseudo_prepid,
                                                   self.sample_name,self.sample_type)
                if production_vcf == None:
                    ## Note this should not happen since the task is only
                    ## yielded upstream if the return is valid
                    raise Exception("Production vcf missing,bug in dependency"
                          "resolution cannot run genotype concordance!!")

                truth_vcf,eval_vcf = self.subset_vcf([production_vcf,self.annotated_vcf_gz])
                if os.path.exists(truth_vcf+'.idx'):
                    os.system("rm {0}.idx".format(truth_vcf))
                if os.path.exists(eval_vcf+'.idx'):
                    os.system("rm {0}.idx".format(eval_vcf))
                with open(self.script,'w') as OUT:
                    OUT.write(self.concordance_cmd.format(config().java,config().gatk,config().ref,eval_vcf,truth_vcf,self.output_file,self.log_file)+'\n')
                proc = subprocess.Popen(self.concordance_cmd.format(config().java,config().gatk,config().ref,eval_vcf,truth_vcf,self.output_file,self.log_file),shell=True)
                proc.wait()
                if proc.returncode: ## Non zero return code
                    raise subprocess.CalledProcessError(proc.returncode,
                                                            self.concordance_cmd)
                remove_cmd = "rm temp.vcf {0} {1}"
                db = get_connection("seqdb")
                cur = db.cursor()
                update_dragen_status(self.pseudo_prepid,self.sample_name,
                                     self.__class__.__name__,cur)
                cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                    pseudo_prepid=self.pseudo_prepid,
                    pipeline_step_id=self.pipeline_step_id))
                db.commit()
                db.close()
            except MySQLdb.Error, e:
                if i == tries - 1:
                    raise Exception("ERROR %d IN CONNECTION: %s" % (e.args[0], e.args[1]))
                else:
                    time.sleep(60) ## Wait a minute before trying again
                    continue
            finally:
                if db.open:
                    db.close()
                break ## Exit from the loop since we have completed the task 
    def subset_vcf(self,vcfs):
        """ Get only PASS variants based on the FILTER field

        Args : vcfs ; List; A list of vcfs to subset, currently only handle 2 vcfs inside the function though
               Assume it to be in the order [truth_vcf,eval_vcf]

        Returns : List; Paths to the output vcf
        """

        subset_cmd = "zcat {0} | awk -v sample={1} -f {2} > {3}"
        out_vcf = [os.path.join(self.sample_dir,'truth.vcf'),os.path.join(self.sample_dir,'eval.vcf')]
        i=0
        for vcf in vcfs:
            cmd = subset_cmd.format(vcf,self.sample_name,config().subset_vcf_awk,out_vcf[i])
            proc = subprocess.Popen(cmd,shell=True)
            proc.wait()
            if proc.returncode: ## Non zero return code
                raise subprocess.CalledProcessError(proc.returncode,cmd)
            i+=1

        return out_vcf

class ContaminationCheck(SGEJobTask):
    """ Run VerifyBamID to check for sample contamination """

    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self,*args,**kwargs):
        super(ContaminationCheck,self).__init__(*args,**kwargs)
        self.seqtype_dir = os.path.join(self.scratch,self.sample_type.upper())
        self.sample_dir = os.path.join(self.seqtype_dir,self.sample_name+'.'+self.pseudo_prepid)
        kwargs["sample_name"] = self.sample_name
        super(ContaminationCheck, self).__init__(*args, **kwargs)
        self.output_stem = os.path.join(self.sample_dir,
                                        "{0}.contamination.raw"
                                        .format(self.sample_name))
        self.output_file = os.path.join(self.sample_dir,
                                        "{0}.{1}.contamination.selfSM"
                                        .format(self.sample_name,self.pseudo_prepid))
        self.annotated_vcf_gz = os.path.join(self.sample_dir,"{0}.{1}.analysisReady.annotated.vcf.gz".format(self.sample_name,self.pseudo_prepid))
        self.snp_vcf = os.path.join(self.sample_dir,"{0}.{1}.snpcontam.vcf.gz".format(self.sample_name,self.pseudo_prepid))
        self.log_dir = os.path.join(self.sample_dir,'logs')
        self.log_file = os.path.join(self.log_dir,
                                     self.sample_name+
                                     '.{0}.samplecontamination.log'.format(self.pseudo_prepid))
        self.script = os.path.join(self.sample_dir,"scripts/{0}.{1}.concordance.metrics.sh".format(self.sample_name,self.pseudo_prepid))
        self.recal_bam = os.path.join(self.sample_dir,self.sample_name+'.{0}.realn.recal.bam'.format(self.pseudo_prepid))
        self.cmd = "{0} --vcf {1} --bam {2} --out {3} --verbose --ignoreRG --maxDepth 1000 --precise >> {4}".format(
            config().verifybamid,config().contam1000g_vcf,self.recal_bam,self.output_stem,self.log_file)
        self.parser_cmd = "awk -f {0} {1}.selfSM > {2}".format(config().transpose_awk,self.output_stem,self.output_file)
        self.remove_cmd = "rm {0}".format(self.snp_vcf)
        self.pipeline_step_id = get_pipeline_step_id(
            self.__class__.__name__,"seqdb")

    def requires(self):
        """
        The requirement for this task is the presence of the analysis ready
        vcf file from the GATK pipeline
        """

        return self.clone(AnnotateVCF)

    def output(self):
        """ The result from this task is the creation of the metrics file """
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

    def work(self):
        with open(self.script,'w') as OUT:
            OUT.write(self.cmd+'\n')
            OUT.write(self.parser_cmd+'\n')

        run_shellcmd("seqdb",self.pipeline_step_id,self.pseudo_prepid,
                     [self.cmd,self.parser_cmd],
                     self.sample_name,self.__class__.__name__)

def update_alignseqfile(location,chgvid,pseudo_prepid):
    """ Update alignseqfile loc in seqdb

    location : str ; the path to the alignseqfileloc
    chgvid : str ; the chgvid for the sample
    pseudo_prepid : str ; the pseudo_prepid for the sample
    """

    update_statement = (""" UPDATE dragen_qc_metrics SET """
                        """ AlignSeqFileLoc = '{0}' WHERE """
                        """ pseudo_prepid = {1} AND """
                        """ CHGVID = '{2}' """.format(
                            location,pseudo_prepid,chgvid))
    try:
        db = get_connection("seqdb")
        cur = db.cursor()
        cur.execute(update_statement)
        db.commit()
        db.close()
    except MySQLdb.Error, e:
        raise Exception("ERROR %d IN CONNECTION: %s" % (e.args[0], e.args[1]))
    finally:
        if db.open:
            db.close()

class UpdateSeqdbMetrics(SGEJobTask):
    """ Populate database with output files """
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    capture_kit_bed = luigi.Parameter()
    sample_type = luigi.Parameter()
    scratch = luigi.Parameter()
    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/tmp/luigi_test"

    def __init__(self,*args,**kwargs):
        super(UpdateSeqdbMetrics,self).__init__(*args,**kwargs)
        self.seqtype_dir = os.path.join(self.scratch,self.sample_type.upper())
        self.sample_dir = os.path.join(self.seqtype_dir,self.sample_name+'.'+self.pseudo_prepid)
        kwargs["sample_name"] = self.sample_name
        super(UpdateSeqdbMetrics, self).__init__(*args, **kwargs)

        self.log_dir = os.path.join(self.sample_dir,'logs')
        self.log_file = os.path.join(self.log_dir,
                                     self.sample_name+
                                     '.{0}.updateseqdb.log'.format(self.pseudo_prepid))
        self.annotated_vcf_gz = os.path.join(self.sample_dir,"{0}.{1}.analysisReady.annotated.vcf.gz".format(self.sample_name,self.pseudo_prepid))
        self.recal_bam = os.path.join(self.sample_dir,"{0}.{1}.realn.recal.bam".format(self.sample_name,self.pseudo_prepid))
        ## Generic query to be used for updates 
        self.update_statement = """ UPDATE {0} SET {1} = '{2}' WHERE CHGVID = '{3}' AND seqType = '{4}' AND pseudo_prepid = '{5}'"""
        self.query_statement = """ SELECT {0} FROM {1} WHERE CHGVID = '{2}' AND seqType = '{3}' AND prepid = {4}"""

        self.issue_contamination_warning = False
        ## The output files from the tasks run 
        self.alignment_out = os.path.join(self.sample_dir,"{0}.{1}.alignment.metrics.txt".format(self.sample_name,self.pseudo_prepid))
        self.cvg_out = os.path.join(self.sample_dir,"{0}.{1}.cvg.metrics.txt".format(self.sample_name,self.pseudo_prepid))
        self.cvg_ccds_out = os.path.join(self.sample_dir,"{0}.{1}.cvg.metrics.ccds.txt".format(self.sample_name,self.pseudo_prepid))
        self.cvg_X_out = os.path.join(self.sample_dir,"{0}.{1}.cvg.metrics.X.txt".format(self.sample_name,self.pseudo_prepid))
        self.cvg_Y_out = os.path.join(self.sample_dir,"{0}.{1}.cvg.metrics.Y.txt".format(self.sample_name,self.pseudo_prepid))
        self.cvg_cs_out = os.path.join(self.sample_dir,"{0}.{1}.cvg.metrics.cs.txt".format(self.sample_name,self.pseudo_prepid))
        self.dup_out = os.path.join(self.sample_dir,"{0}.{1}.duplicates.txt".format(self.sample_name,self.pseudo_prepid))
        self.variant_call_out = os.path.join(self.sample_dir,"{0}.{1}.variant_calling_summary_metrics".format(self.sample_name,self.pseudo_prepid))
        self.geno_concordance_out = os.path.join(self.sample_dir,"{0}.{1}.genotype_concordance_metrics".format(self.sample_name,self.pseudo_prepid))
        self.contamination_out = os.path.join(self.sample_dir,"{0}.{1}.contamination.selfSM".format(self.sample_name,self.pseudo_prepid))
        ## The qc table to update
        self.qc_table = "dragen_qc_metrics"
        ## The new lean equivalent of seqdbClone 
        self.master_table = "seqdbClone"
        self.DBID,self.prepID = getDBIDMaxPrepID(self.pseudo_prepid)
        ## Get the id for this pipeline step
        self.pipeline_step_id = get_pipeline_step_id(
            self.__class__.__name__,"seqdb")

        self.raw_vcfs = os.path.join(self.sample_dir,"{0}.{1}.raw*vcf*".format(self.sample_name,self.pseudo_prepid))
        self.indel_vcfs = os.path.join(self.sample_dir,"{0}.{1}.indel*".format(self.sample_name,self.pseudo_prepid))
        self.snp_vcfs = os.path.join(self.sample_dir,"{0}.{1}.snp*".format(self.sample_name,self.pseudo_prepid))
        self.analysis_vcfs = os.path.join(self.sample_dir,"{0}.{1}.analysisReady.vcf*".format(self.sample_name,self.pseudo_prepid))
        self.fixed_vcf = os.path.join(self.sample_dir,"{0}.{1}.analysisReady.fixed.vcf".format(self.sample_name,self.pseudo_prepid))
        self.tmp_vcf = os.path.join(self.sample_dir,"{0}.{1}.tmp.vcf".format(self.sample_name,self.pseudo_prepid))
        self.files_to_remove = [self.raw_vcfs,self.indel_vcfs,self.snp_vcfs,self.analysis_vcfs,self.fixed_vcf,self.tmp_vcf]

    def requires(self):
        """
        The requirement for this task
        the outputfile from different tasks should
        be present
        """

        yield self.clone(AnnotateVCF)
        yield self.clone(GQBinning)
        yield self.clone(CvgBinning)
        yield self.clone(AlignmentMetrics)
        yield self.clone(RunCvgMetrics)
        yield self.clone(DuplicateMetrics)
        yield self.clone(VariantCallingMetrics)
        yield self.clone(ContaminationCheck)
        if get_productionvcf(self.pseudo_prepid,self.sample_name,self.sample_type) != None:
            yield self.clone(GenotypeConcordance)

    def output(self):
        return SQLTarget(pseudo_prepid=self.pseudo_prepid,
                         pipeline_step_id=self.pipeline_step_id)

    def work(self):
        tries = 5
        for i in range(tries):
            try:
                ## Get db connection and cursor
                db = get_connection("seqdb")
                cur = db.cursor()
                ## Update pipeline start time
                cur.execute(UPDATE_PIPELINE_STEP_SUBMIT_TIME.format(
                    pseudo_prepid=self.pseudo_prepid,
                    pipeline_step_id=self.pipeline_step_id))
                db.commit()
                db.close()

                self.LOG_FILE = open(self.log_file,'w')
                self.add_sample_to_qc()
                self.update_alignment_metrics()
                self.update_coverage_metrics('all')
                self.update_coverage_metrics('ccds')
                self.update_coverage_metrics('X')
                self.update_coverage_metrics('Y')
                if self.sample_type.upper() != 'GENOME':
                    self.update_coverage_metrics('cs')
                self.update_duplicates()
                self.update_variantcalling_metrics()
                if get_productionvcf(self.pseudo_prepid,self.sample_name,self.sample_type) != None:
                    self.update_genotype_concordance_metrics()
                self.update_contamination_metrics()
                self.update_seqgender()
                self.update_qc_message()
                ## Update alignseqfile location here to keep track of scratch location
                update_alignseqfile(self.seqtype_dir,self.sample_name,self.pseudo_prepid)
                ## Remove all temp files
                for tmp_file in self.files_to_remove:
                    os.system("rm {0}".format(tmp_file))

                db = get_connection("seqdb")
                cur = db.cursor()
                ## Update dragen_status
                update_dragen_status(self.pseudo_prepid,self.sample_name,
                                     self.__class__.__name__,cur)
                ## Update pipeline finish time
                cur.execute(UPDATE_PIPELINE_STEP_FINISH_TIME.format(
                    pseudo_prepid=self.pseudo_prepid,
                    pipeline_step_id=self.pipeline_step_id))
                db.commit()
                db.close()
                break
            except MySQLdb.Error, e:
                if i > tries - 1:
                    raise Exception("ERROR %d IN CONNECTION: %s" %(e.args[0],e.args[1]))
                else:
                    time.sleep(60)
                    continue
            finally:
                self.LOG_FILE.close()
                if db.open:
                    db.close()

    def add_sample_to_qc(self):
        """
        Initialize the sample in the qc table, use update statements later to update qc metrics

        Args : Nothing

        Returns : Nothing
        """

        query = """INSERT INTO {0} (CHGVID,seqType,pseudo_prepid) VALUES ('{1}','{2}',{3})""".format(self.qc_table,self.sample_name,self.sample_type,self.pseudo_prepid)
        try:
            db = get_connection("seqdb")
            cur = db.cursor()
            cur.execute(query)
            db.commit()
            db.close()
        except MySQLdb.IntegrityError as e:
            if not e[0] == 1062: ## The entry is present in the qc table
                raise Exception("ERROR %d IN CONNECTION: %s" %(e.args[0],e.args[1]))
            else:
                print >> self.LOG_FILE,"This chgvid,prepid,seqtype was already present in the qc table, this occurs upon retry of the task, not a bug"

    def check_qc(self):
        """
        Check if all the metrics meet the thresholds,will call individual
        functions for checking the appropriate metric

        Args : Nothing

        Returns : Boolean ; True/False
        """

        if (self.check_alignment() and self.check_duplicates() and self.check_variantcalling() and self.check_coverage() and self.check_contamination()):
            return True
        else:
            return False

    def update_qc_message(self):
        """ Update qc message in statusT based on qc check, will call the check_qc
        function to perform the check """

        self.failed_qc = "QC review needed"
        self.passed_qc = "Passed Bioinfo QC"
        message=[]
        if self.check_qc():
            message.append(self.passed_qc)
        else: ## Check individual qc checks again for updating qc message
            message.append(self.failed_qc)
            if self.issue_contamination_warning == True: ## Check for contamination warning
                message.append("Warning sample contamination is high, but below qc fail threshold")
            if not self.check_alignment():
                message.append("Failed Alignment Check")
            if not self.check_duplicates():
                message.append("Failed Duplicate Check")
            if not self.check_variantcalling():
                message.append("Failed VariantCallingCheck")
            if not self.check_coverage():
                message.append("Failed Coverage Check")
            if not self.check_contamination():
                message.append("Failed Contamination Check")
            #if not self.check_gender():
                #message.append("Failed Gender Check")

        final_message = ';'.join(message)
        self.updatedatabase(self.qc_table,'QCMessage',final_message)


    def updatedatabase(self,table_name,db_field,val):
        """
        Function to update qc fields in seqdbClone
        Will raise an exception if there was a connection
        error

        Args : table_name ; String ; The database table to update
               db_field ; String ; The database field to update
               val ; String ; The value to update with

        Return : Does not return anything
        """

        tries = 1
        for i in range(tries):
            try:
                db = get_connection("seqdb")
                cur = db.cursor()
                cur.execute(self.update_statement.format(table_name,db_field,val,self.sample_name,self.sample_type,self.pseudo_prepid))
                db.commit()
                db.close()
            except MySQLdb.Error, e:
                if i == tries - 1:
                    raise Exception("ERROR %d IN CONNECTION: %s" % (e.args[0], e.args[1]))
                else:
                    time.sleep(60) ## Wait a minute before trying again
                    continue
            finally:
                if db.open:
                    db.close()
                break ## Exit the loop since the update is complete

    def get_metrics(self,query):
        """
        Queries seqdbClone and returns the result, will always return a tuple
        since we do not know apriori how many how many fields are being needed

        Args : query ; String ; The query to be executed
        Returns : tuple ; The result from the query
        """

        tries = 1
        for i in range(tries):
            try:
                db = get_connection("seqdb")
                cur = db.cursor()
                cur.execute(query)
                db_val = cur.fetchall()
                db.close()
                return db_val ## Pass control back
            except MySQLdb.Error, e:
                if i == tries - 1:
                    raise Exception("ERROR %d IN CONNECTION: %s" % (e.args[0], e.args[1]))
                else:
                    time.sleep(60) ## Wait a minute before trying again
                    continue

    def update_alignment_metrics(self):
        """
        Function to update the Alignment Metrics
        Makes a call to the generic updatedatabase function with the
        appropriate database field and value to update with

        Args : None
        Returns : Does not return anything
        """

        self.alignment_parse = {'TOTAL_READS':qcmetrics().total_reads,'PCT_PF_READS_ALIGNED':qcmetrics().pct_reads_aligned,'PF_MISMATCH_RATE':qcmetrics().pct_mismatch_rate}

        with open(self.alignment_out) as OUTFILE:
            for line in OUTFILE:
                contents = line.strip().strip('\n').rstrip(' ').split(' ')
                if len(contents) == 4:
                    field = contents[0]
                    val = contents[-1]
                    if field in self.alignment_parse.keys():
                        db_field = self.alignment_parse[field]
                        self.updatedatabase(self.qc_table,db_field,val)

    def update_coverage_metrics(self,file_type):
        """
        Function to update the Coverage Metrics
        Makes a call to the generic updatedatabase function with the
        appropriate database field and value to update with

        Args : file_type ; String ; Can be one of the following : ['all','ccds','X','Y','cs']
        Returns : Does not return anything
        """

        self.cvg_exome_parse_ccds = {'mean':qcmetrics().mean_ccds_cvg,'%_bases_above_5':qcmetrics().pct_ccds_bases5X,'%_bases_above_10':qcmetrics().pct_ccds_bases10X,'%_bases_above_15':qcmetrics().pct_ccds_bases15X,'%_bases_above_20':qcmetrics().pct_ccds_bases20X}
        self.cvg_exome_parse = {'mean':qcmetrics().mean_cvg,'granular_median':qcmetrics().median_cvg,'%_bases_above_5':qcmetrics().pct_bases5X,'%_bases_above_10':qcmetrics().pct_bases10X,'%_bases_above_15':qcmetrics().pct_bases15X,'%_bases_above_20':qcmetrics().pct_bases20X}
        self.cvg_exome_parse_X = {'mean':qcmetrics().mean_X_cvg}
        self.cvg_exome_parse_Y = {'mean':qcmetrics().mean_Y_cvg}
        self.cvg_exome_parse_cs = {'ON_BAIT_VS_SELECTED':qcmetrics().capture_specificity}
        self.cvg_wgs_parse_ccds = {'MEAN_COVERAGE':qcmetrics().mean_cvg,'PCT_5X':qcmetrics().pct_bases5X,'PCT_10X':qcmetrics().pct_bases10X,'PCT_15X':qcmetrics().pct_bases15X,'PCT_20X':qcmetrics().pct_bases20X}
        self.cvg_wgs_parse_X = {'MEAN_COVERAGE':qcmetrics().mean_X_cvg}
        self.cvg_wgs_parse_Y = {'MEAN_COVERAGE':qcmetrics().mean_Y_cvg}
        self.cvg_wgs_parse = {'MEAN_COVERAGE':qcmetrics().mean_cvg,'MEDIAN_COVERAGE':qcmetrics().median_cvg,'PCT_5X':qcmetrics().pct_bases5X,'PCT_10X':qcmetrics().pct_bases10X,'PCT_15X':qcmetrics().pct_bases15X,'PCT_20X':qcmetrics().pct_bases20X}

        if self.sample_type.upper() == 'GENOME':
            if file_type == 'ccds':
                metrics_hash = self.cvg_wgs_parse_ccds
                metrics_file = self.cvg_ccds_out
            elif file_type == 'X':
                metrics_hash = self.cvg_wgs_parse_X
                metrics_file = self.cvg_X_out
            elif file_type == 'Y':
                metrics_hash = self.cvg_wgs_parse_Y
                metrics_file = self.cvg_Y_out
            elif file_type == 'all':
                metrics_hash = self.cvg_wgs_parse
                metrics_file = self.cvg_out
        else: ## Exome or custom capture
            if file_type == 'ccds':
                metrics_hash = self.cvg_exome_parse_ccds
                metrics_file = self.cvg_ccds_out
            elif file_type == 'X':
                metrics_hash = self.cvg_exome_parse_X
                metrics_file = self.cvg_X_out
            elif file_type == 'Y':
                metrics_hash = self.cvg_exome_parse_Y
                metrics_file = self.cvg_Y_out
            elif file_type == 'cs':
                metrics_hash = self.cvg_exome_parse_cs
                metrics_file = self.cvg_cs_out
            elif file_type == 'all':
                metrics_hash = self.cvg_exome_parse
                metrics_file = self.cvg_out

        with open(metrics_file) as OUTFILE:
            for line in OUTFILE:
                contents = line.strip().strip('\n').rstrip(' ').split(' ')
                if len(contents) > 1:
                    field = contents[0]
                    val = contents[1]
                    if field in metrics_hash.keys():
                        db_field = metrics_hash[field]
                        self.updatedatabase(self.qc_table,db_field,val)
                        if file_type == 'all':
                            if db_field == 'OverallCoverage':
                                mean_cvg = float(val)
                            if db_field == 'MedianCoverage':
                                median_cvg = float(val)

            if file_type == 'all':
                db_field = qcmetrics().mean_median_cvg
                val = mean_cvg/median_cvg
                self.updatedatabase(self.qc_table,db_field,val)


    def update_duplicates(self):
        """
        Function to update the Variant Calling Metrics
        Makes a call to the generic updatedatabase function with the
        appropriate database field and value to update with

        Args : None
        Returns : Does not return anything
        """

        with open(self.dup_out,'r') as OUTFILE:
            for line in OUTFILE:
                line = line.strip('\n')
                self.updatedatabase(self.qc_table,qcmetrics().pct_duplicate_reads,line)

    def update_variantcalling_metrics(self):
        """
        Function to update the Variant Calling Metrics
        Makes a call to the generic updatedatabase function with the
        appropriate database field and value to update with

        Args : None
        Returns : Does not return anything
        """

        ## Need to calculate overall titv since that is not available in the output file
        self.variant_call_parse = {'TOTAL_SNPS':qcmetrics().total_snps,'PCT_DBSNP':qcmetrics().pct_dbsnp_snps,'TOTAL_INDELS':qcmetrics().total_indels,'PCT_DBSNP_INDELS':qcmetrics().pct_dbsnp_indels}
        temp_hash = {}
        flag = 0 
        with open(self.variant_call_out,'r') as OUTFILE:
            for line in OUTFILE:
                field,val = line.strip().strip('\n').rstrip(' ').split(' ')[0:2]
                if field in self.variant_call_parse.keys():
                    db_field = self.variant_call_parse[field]
                    self.updatedatabase(self.qc_table,db_field,val)
                temp_hash[field] = val

            titv = (((int(temp_hash['NOVEL_SNPS'])*float(temp_hash['NOVEL_TITV']) + (int(temp_hash['NUM_IN_DB_SNP'])*float(temp_hash['DBSNP_TITV']))))/int(temp_hash['TOTAL_SNPS']))
            self.updatedatabase(self.qc_table,qcmetrics().titv,titv)
            ## Note i make calls to helper functions in the same class to calculate the hom and het counts for different conditions 
            homhet_ratio = float(self.get_hom_count())/float(self.get_het_count())
            self.updatedatabase(self.qc_table,qcmetrics().homhet_ratio,homhet_ratio)
            snv_homhet_ratio = float(self.get_hom_count(False,False,True))/float(self.get_het_count(False,False,True))
            self.updatedatabase(self.qc_table,qcmetrics().snv_homhet_ratio,snv_homhet_ratio)
            indel_homhet_ratio = float(self.get_hom_count(False,True,False))/float(self.get_het_count(False,True,False))
            self.updatedatabase(self.qc_table,qcmetrics().indel_homhet_ratio,indel_homhet_ratio)
            x_homhet_ratio = float(self.get_hom_count(True))/float(self.get_het_count(True))
            self.updatedatabase(self.qc_table,qcmetrics().x_homhet_ratio,x_homhet_ratio)
    def update_genotype_concordance_metrics(self):
        """
        Function to update the Genotype Concordance Metrics
        Makes a call to the generic updatedatabase function with the
        appropriate database field and value to update with

        Args : None
        Returns : Does not return anything
        """

        self.geno_concordance_parse = {'ALLELES_MATCH':1,'ALLELES_DO_NOT_MATCH':1,'EVAL_ONLY':1,'TRUTH_ONLY':1}
        self.temp_genoconcordance = os.path.join(self.sample_dir,"temp_geno_concordance.txt")
        cmd = """cat {0} | grep -A 2 "#:GATKTable:SiteConcordance_Summary:Site-level summary statistics"| grep -v "#" | awk -f {1} > {2}""".format(self.geno_concordance_out,config().transpose_awk,self.temp_genoconcordance)
        proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        proc.wait()
        if proc.returncode:
            print >> self.LOG_FILE,subprocess.CalledProcessError(proc.returncode,self.get_het_cmd)
        with open(self.temp_genoconcordance,'r') as OUTFILE:
            for line in OUTFILE:
                field,val = line.strip().strip('\n').rstrip(' ').split(' ')[0:2]
                if field in self.geno_concordance_parse.keys():
                    self.geno_concordance_parse[field] = int(val)

            concordance_metric = (float(self.geno_concordance_parse['ALLELES_MATCH']))/(self.geno_concordance_parse['ALLELES_DO_NOT_MATCH'] + self.geno_concordance_parse['ALLELES_MATCH'] + self.geno_concordance_parse['EVAL_ONLY'] + self.geno_concordance_parse['TRUTH_ONLY'])
            self.updatedatabase(self.qc_table,qcmetrics().concordance,concordance_metric)

    def update_contamination_metrics(self):
        """
        Function to update the Contamination Metrics
        Makes a call to the generic updatedatabase function with the
        appropriate database field and value to update with

        Args : None
        Returns : Does not return anything
        """

        with open(self.contamination_out) as OUTFILE:
            for line in OUTFILE:
                field,val = line.strip().strip('\n').rstrip(' ').split(' ')[0:2]
                if field == 'FREEMIX':
                    db_field = qcmetrics().contamination_value
                    self.updatedatabase(self.qc_table,db_field,val)
                    ## No need to loop further
                    break

    def update_seqgender(self):
        """
        Calculates the X/Y coverage for Exomes and Genomes and updates the seq gender
        based on existing rules : https://redmine.igm.cumc.columbia.edu/projects/biopipeline/wiki/Gender_checks
        For custom capture regions the rules are more complicated. 
        Also updates the database with het/hom ratio on X chromosome, the num_homX and the num_hetX
        Args : Nothing
        Returns : String; One of the following values : [Male,Female,Ambiguous]
        """

        if self.sample_type.upper() in ['EXOME','GENOME']:
            query = """SELECT {0} FROM {1} WHERE CHGVID = '{2}' AND seqType = '{3}' AND pseudo_prepid = '{4}'"""
            result = self.get_metrics(query.format(qcmetrics().mean_X_cvg,self.qc_table,self.sample_name,self.sample_type,self.pseudo_prepid))
            if len(result) == 0:
                seq_gender = 'Ambiguous'
            else:
                X_cvg = float(result[0][0])
                result = self.get_metrics(query.format(qcmetrics().mean_Y_cvg,self.qc_table,self.sample_name,self.sample_type,self.pseudo_prepid))
                if len(result) == 0:
                    seq_gender = 'Ambiguous'
                else:
                    Y_cvg = float(result[0][0])
                    if X_cvg/Y_cvg < 2:
                        seq_gender = 'M'
                    elif X_cvg/Y_cvg > 5:
                        seq_gender = 'F'
                    else:
                        seq_gender = 'Ambiguous'

        else : ## Custom Capture Samples
            het = self.get_het_count(sex=True)
            hom = self.get_hom_count(sex=True)
            if het == 0:
                seq_gender = 'M'
            elif hom == 0:
                seq_gender = 'F'
            elif float(het)/hom < 0.26:
                seq_gender = 'M'
            elif float(het)/hom > 0.58:
                seq_gender = 'F'
            else:
                seq_gender = 'Ambiguous'

        ## Update the database
        if Y_cvg != 0:
            x_y_cvg = X_cvg/Y_cvg
            self.updatedatabase(self.qc_table,qcmetrics().xy_avg_ratio,x_y_cvg)
        self.updatedatabase(self.qc_table,'SeqGender',seq_gender)

    def check_alignment(self):
        """
        Checks if the alignment metrics matches the appropriate thresholds
        Currently using only the perc_reads_aligned, will use the get_metrics
        function to get the requisitve metric(s) required.

        Args : Nothing
        Returns : Bool; True/False
        """

        query = """SELECT {0} FROM {1} WHERE CHGVID = '{2}' AND seqType = '{3}' AND pseudo_prepid = '{4}'""".format(qcmetrics().pct_reads_aligned,self.qc_table,self.sample_name,self.sample_type,self.pseudo_prepid)
        result = self.get_metrics(query)
        perc_reads_aligned = float(result[0][0])
        return (perc_reads_aligned >= 0.70)

    def check_coverage(self):
        """
        Checks if the coverage metrics matches the appropriate thresholds
        Currently not using anything

        Args : Nothing
        Returns : Bool; True/False
        """

        query = """SELECT {0},{1} FROM {2} WHERE CHGVID = '{3}' AND seqType = '{4}' AND pseudo_prepid = '{5}'""".format(qcmetrics().pct_bases5X,qcmetrics().mean_cvg,self.qc_table,self.sample_name,self.sample_type,self.pseudo_prepid)
        result = self.get_metrics(query)
        pct_bases_5X = float(result[0][0])
        mean_cvg = float(result[0][1])
        return (pct_bases_5X >= 0.9 and mean_cvg >= 25)

    def check_duplicates(self):
        """
        Checks if the duplicate metrics matches the appropriate thresholds
        Currently using only the perc_duplicate_reads, will use the get_metrics
        function to get the requisitve metric(s) required.

        Args : Nothing
        Returns : Bool; True/False
        """

        query = """SELECT {0} FROM {1} WHERE CHGVID = '{2}' AND seqType = '{3}' AND pseudo_prepid = '{4}'""".format(qcmetrics().pct_duplicate_reads,self.qc_table,self.sample_name,self.sample_type,self.pseudo_prepid)
        result = self.get_metrics(query)
        perc_duplicate_reads = float(result[0][0])
        if self.sample_type.upper() == 'GENOME':
            return (perc_duplicate_reads <= 20)

        else:
            return (perc_duplicate_reads <= 30)


    def check_variantcalling(self):
        """
        Checks if the concordance metrics matches the appropriate thresholds
        Currently using only the , will use the get_metrics
        function to get the requisitve metric(s) required.

        Args : num_snvs ; type:can handle str and int ; the number of snvs in the vcf
        Returns : Bool; True/False
        """

        query = """SELECT {0} FROM {1} WHERE CHGVID = '{2}'AND seqType = '{3}' AND pseudo_prepid = '{4}'""".format(qcmetrics().total_snps,self.qc_table,self.sample_name,self.sample_type,self.pseudo_prepid)
        result = self.get_metrics(query)
        num_snvs = int(result[0][0])

        if self.sample_type.upper() == 'GENOME':
            return (int(num_snvs) > 3000000)
        else:
            return (int(num_snvs) > 100000)

    def check_concordance(self):
        """
        Checks if the concordance metrics matches the appropriate thresholds
        Currently not using anything

        Args : Nothing
        Returns : Bool; True/False
        """

        return True

    def check_contamination(self):
        """
        Checks if the contamination metrics matches the appropriate thresholds
        Currently using only the , will use the get_metrics
        function to get the requisitve metric(s) required.

        Args : Nothing
        Returns : Bool; True/False
        """

        query = """SELECT {0} FROM {1} WHERE CHGVID = '{2}' AND seqType = '{3}' AND pseudo_prepid = '{4}'""".format(qcmetrics().contamination_value,self.qc_table,self.sample_name,self.sample_type,self.pseudo_prepid)
        result = self.get_metrics(query)
        contamination_value = float(result[0][0])
        if contamination_value > 0.05 and contamination_value < 0.08:
            self.issue_contamination_warning = True
            return True
        elif contamination_value > 0.08:
            return False
        else:
            return True

    def check_gender(self):
        """
        Query seqdbClone and check whether the seq and selfdecl gender matchup
        Args : Nothing
        Returns : Bool; True/False
        """

        query = """SELECT {0} FROM {1} WHERE CHGVID = '{2}' AND seqType = '{3}' AND pseudo_prepid = '{4}'"""
        q1 = query.format('SelfDeclGender',self.master_table,self.sample_name,self.sample_type,self.pseudo_prepid)
        result = self.get_metrics(q1)
        selfdecl_gender = result[0][0]
        q2 = query.format('SeqGender',self.master_table,self.sample_name,self.sample_type,self.pseudo_prepid)
        result = self.get_metrics(q2)
        seq_gender = result[0][0]
        return (selfdecl_gender == seq_gender)

    def get_hom_count(self,sex=False,indel=True,snv=True,**kwargs):
        """
        Returns number of homozygous sites in the vcf
        matching certain conditions
        Try to use JEXL in the future: https://software.broadinstitute.org/gatk/guide/article?id=1255

        Args : sex ; Bool ; Calculate on the X chromosome, by default False
               indel ; Bool ; If all the full vcf , then whether to restrict only to indels
               snv ; Bool ; If the full vcf, then whether to restrict only to snv
               Note : only one of snv or indel can be False

        Returns : Int; Count for the number of hom sites
        """

        ## Override default arguments if they are present 
        if 'sex' in kwargs:
            sex = kwargs['sex']
        if 'indel' in kwargs:
            indel = kwargs['indel']
        if 'snv' in kwargs:
            snv = kwargs['snv']

        if sex == True: ## We will always use both snps and indels on the X chromosome
            self.get_hom_cmd = """ %s -f "MQ > 39 & QD > 1" -g "GQ > 19" %s | grep -v "#" | grep "^X" | awk '{print $NF}'|awk -F ":" '{print $1}' | awk -F "/" '$1==$2{print}'| wc -l"""%(config().vcffilter,self.annotated_vcf_gz)
        elif indel == False: ## Only SNPs
            self.get_hom_cmd = """ zcat %s | grep -v "#" | grep "PASS" | awk '{if (length($4) == length($5)){print $NF}}'|awk -F ":" '{print $1}' | awk -F "/" '$1==$2{print}'| wc -l"""%(self.annotated_vcf_gz)
        elif snv == False: ## Only Indels
            self.get_hom_cmd = """ zcat %s | grep -v "#" | grep "PASS" | awk '{if (length($4)!=length($5)){print $NF}}' |awk -F ":" '{print $1}' | awk -F "/" '$1==$2{print}'| wc -l"""%(self.annotated_vcf_gz)
        else: ## Both SNPs and Indels
            self.get_hom_cmd = """ zcat %s | grep -v "#" | grep "PASS" | awk '{print $NF}'|awk -F ":" '{print $1}' | awk -F "/" '$1==$2{print}'| wc -l"""%(self.annotated_vcf_gz)

        proc = subprocess.Popen(self.get_hom_cmd,shell=True,stdout=subprocess.PIPE)
        proc.wait()
        if proc.returncode:
            print >> self.LOG_FILE,subprocess.CalledProcessError(proc.returncode,self.get_het_cmd)

        return int(proc.stdout.read().strip('\n'))

    def get_het_count(self,sex=False,indel=True,snv=True):
        """
        Returns number of heterozygous sites in the vcf
        matching certain conditions
        Try to use JEXL in the future : https://software.broadinstitute.org/gatk/guide/article?id=1255

        Args : sex ; Bool ; Calculate on the X chromosome, by default False
               indel ; Bool ; If all the full vcf , then whether to restrict only to indels
               snv ; Bool ; If the full vcf, then whether to restrict only to snv
               Note : only one of snv or indel can be False

        Returns : Int; Count for the number of het sites
        """

        if sex == True: ## We will always use both snps and indels on the X chromosome, there are some additional filtering criteria which are being applied, see here for more details : https://redmine.igm.cumc.columbia.edu/projects/biopipeline/wiki/Gender_checks
            self.get_het_cmd = """ %s -f "MQ > 39 & QD > 1" -g "GQ > 19" %s| grep -v "#" | grep "^X" | awk '{print $NF}'|awk -F ":" '{print $1}' | awk -F "/" '$1!=$2{print}'| wc -l"""%(config().vcffilter,self.annotated_vcf_gz)
        elif indel == False: ## Only SNPs
            self.get_het_cmd = """ zcat %s | grep -v "#" | grep "PASS" | awk '{if (length($4) == length($5)){print $NF}}'|awk -F ":" '{print $1}' | awk -F "/" '$1!=$2{print}'| wc -l"""%(self.annotated_vcf_gz)
        elif snv == False: ## Only Indels
            self.get_het_cmd = """ zcat %s | grep -v "#" | grep "PASS" | awk '{if (length($4) != length($5)){print $NF}}' | awk -F ":" '{print $1}' | awk -F "/" '$1!=$2{print}'| wc -l"""%(self.annotated_vcf_gz)
        else: ## Both SNPs and Indels
            self.get_het_cmd = """ zcat %s | grep -v "#" | grep "PASS" | awk '{print $NF}'|awk -F ":" '{print $1}' | awk -F "/" '$1!=$2{print}'| wc -l"""%(self.annotated_vcf_gz)

        proc = subprocess.Popen(self.get_het_cmd,shell=True,stdout=subprocess.PIPE)
        proc.wait()
        if proc.returncode:
            print >> self.LOG_FILE,subprocess.CalledProcessError(proc.returncode,self.get_het_cmd)

        return int(proc.stdout.read().strip('\n'))
