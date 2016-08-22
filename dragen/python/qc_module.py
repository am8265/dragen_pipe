#!/nfs/goldstein/software/python2.7.7/bin/python

import os
import shlex
import sys
import subprocess
import luigi
import MySQLdb
from luigi.contrib.sge import SGEJobTask, LocalSGEJobTask
from luigi.util import inherits
import logging
import time
import utils 


class config(luigi.Config):
    """
    config class for instantiating parameters for this pipeline
    the values are read from luigi.cfg in the current folder
    """

    java = luigi.Parameter()
    picard = luigi.Parameter()
    ref = luigi.Parameter()
    seqdict_file = luigi.Parameter()
    bed_file = luigi.Parameter()
    target_file = luigi.Parameter()
    create_targetfile = luigi.BooleanParameter()
    bait_file = luigi.Parameter()
    bait_file_X = luigi.Parameter()
    bait_file_Y = luigi.Parameter()
    python_path = luigi.Parameter()
    relatedness_refs = luigi.Parameter()
    sampleped_loc = luigi.Parameter()
    relatedness_markers = luigi.Parameter()
    bedtools_loc = luigi.Parameter()
    pypy_loc = luigi.Parameter()
    binner_loc  = luigi.Parameter()
    dbsnp = luigi.Parameter()
    cnf_file = luigi.Parameter()
    max_mem = luigi.IntParameter()
    

class db(luigi.Config):
    """
    Database config variable will be read from the
    db section of the config file
    """
    
    cnf = luigi.Parameter()
    seqdb_group = luigi.Parameter()
    dragen_group = luigi.Parameter()

class RootTask(luigi.WrapperTask):
    """
    Wrapper Task
    """

    ## Read in the bam file from the cmd line just to make sure
    ## to ensure uniqueness of the RootTask 
    bam = luigi.Parameter()
    
    
    def requires(self):
        self.output_file = "{0}/{1}.cvg.metrics.txt".format(self.scratch,self.sample_name)
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
    All Picard outputs seem to follow a standard format
    A simple shell command is used for creating a standard processed
    output file to extract values and load into a database
    Information is represented as a two column output,first is the
    field name and the second the value of that field
    """


    picard_output = luigi.Parameter()
    
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"

    def __init__(self,*args,**kwargs):
        """
        Initialize Class Variables
        """

        super(ParsePicardOutputs,self).__init__(*args,**kwargs)
        self.output_file = "{0}.processed".format(self.picard_output)
        

        ## Define the command to be run 
        self.cmd = """cat %s | grep -v '^#'| awk -f /home/rp2801/git/dragen_pipe/dragen/python/transpose.awk > %s """%(self.picard_output,self.output_file)

    def requires(self):
        """ 
        Dependency is the presence of the picard outputfile 
        """
        
        return MyExtTask(self.picard_output)

    def output(self):
        """
        Output from this task
        """

        return luigi.LocalTarget(self.output_file)

    def work(self):
        """
        Just a simple shell command is enough!
        """
        
        os.system(self.cmd)
        


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
        if not os.path.isdir(self.scratch): ## Recursively create the directory if it doesnt exist
            os.makedirs(self.scratch)
            


        ## Define on the fly parameters 
        self.output_file = os.path.join(self.scratch,self.sample_name + ".cvg.metrics.txt")
        self.raw_output_file = os.path.join(self.scratch,self.sample_name + ".cvg.metrics.raw.txt")
        if self.x_y:
            self.output_file_X = os.path.join(self.scratch,self.sample_name+ ".cvg.metrics.X.txt")
            self.raw_output_file_X = os.path.join(self.scratch,self.sample_name + ".cvg.metrics.X.raw.txt")
            self.output_file_Y = os.path.join(self.scratch,self.sample_name+ ".cvg.metrics.Y.txt")
            self.raw_output_file_Y = os.path.join(self.scratch,self.sample_name + ".cvg.metrics.Y.raw.txt")

            
        if self.create_targetfile:
            self.bed_stem = utils.get_filename_from_fullpath(config().bed_file)
            self.target_file = os.path.join(self.scratch,self.bed_stem)

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
    scratch = luigi.Parameter()

    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/rp2801/genome_cvg_temp"


    def __init__(self,*args,**kwargs):
        super(CreatePed,self).__init__(*args,**kwargs)
        if not os.path.isdir(self.scratch): ## Recursively create the directory if it doesnt exist
            os.makedirs(self.scratch)

        self.output_file = self.sample_name+'.'+self.ped_type+'.ped'
        self.output_ped = os.path.join(self.scratch,self.output_file)
                

        if self.ped_type == 'relatedness':
            self.sampletoped_cmd = (
                                    "{0} {1} --scratch {2} --vcffile {3} --bamfile {4}" 
                                    " --markers {5} --refs {6} --sample {7} --fid {8}"
                                    " --seq_type {9} --ped_type {10} --output {11}"
                                    .format(config().python_path,config().sampleped_loc,
                                    self.scratch,self.vcf,self.bam,
                                    config().relatedness_markers,
                                    config().relatedness_refs,self.sample_name,
                                    self.family_id,self.seq_type,self.ped_type,
                                    self.scratch)
                                   )
            
        elif self.ped_type == 'ethnicity':
            self.sampletoped_cmd = (
                                    "{0} {1} --scratch {2} --vcffile {3} --bamfile {4}"
                                    " --markers {5} --refs {6} --sample {7} --fid {8}"
                                    " --seq_type {9} --ped_type {10} --output {11}" 
                                    .format(config().python_path,config().sampleped_loc,
                                    self.scratch,self.vcf,self.bam,
                                    config().ethnicity_markers,
                                    config().ethnicity_refs,self.sample_name,
                                    self.family_id,self.seq_type,self.ped_type,
                                    self.scratch)
                                  )


    def requires(self):
        """
        The dependencies for this task are the existence of the relevant files
        """
        
        return [MyExtTask(self.bam),MyExtTask(self.vcf)]


    def output(self):
        """
        Output from this task is a ped file named with correct type
        """
           
        return luigi.LocalTarget(self.output_ped)


    def work(self):
        """
        To run this task, the sampletoped_fixed.py script is called
        """
        os.system(self.sampletoped_cmd)


       
        
class CreateGenomeBed(SGEJobTask):
    """
    Use bedtools to create the input bed file for the coverage binning script
    """

    ## System Parameters
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/rp2801/genome_cvg_temp/"
    no_tarball=True
    bam = luigi.Parameter()
    sample = luigi.Parameter()
    scratch = luigi.Parameter()
    
    #no_tarball = luigi.Parameter()
    
    def __init__(self,*args,**kwargs):
        super(CreateGenomeBed,self).__init__(*args,**kwargs)

        #self.bam = args[0]
        #self.sample = args[1]
        #self.scratch = args[2]
        #self.poll_time = luigi.Parameter(default=5,description="The time to wait to run qstat")
        
        if not os.path.isdir(self.scratch): ## Recursively create the directory if it doesnt exist
            os.makedirs(self.scratch)

        
        self.genomecov_bed = os.path.join(self.scratch,self.sample+'.genomecvg.bed')

        self.genomecov_cmd = "{0} genomecov -bga -ibam {1} > {2}".format(
                              config().bedtools_loc,self.bam,
                              self.genomecov_bed)

        

        self.gzip_cmd = "bgzip {0}".format(self.genomecov_bed)
        self.gzipped_genomecov_bed = self.genomecov_bed+'.gz'
        self.tabix_cmd = "/nfs/goldstein/software/bin/tabix -p bed {0}".format(self.gzipped_genomecov_bed)
        self.tbi_file = self.gzipped_genomecov_bed+'.tbi'
          
    def requires(self):
        """
        The dependency is the presence of the bam file
        """
        
        return MyExtTask(self.bam)

    def output(self):
        """
        Output is the genomecov bed file
        """

        #return [luigi.LocalTarget(self.genomecov_bed),luigi.LocalTarget(self.gzipped_genomecov_bed),luigi.LocalTarget(self.tbi_file)]
        return [luigi.LocalTarget(self.genomecov_bed)]

                
    def work(self):
        """
        Run the bedtools cmd
        """

        os.system(self.genomecov_cmd)
        #os.system(self.gzip_cmd)
        #os.system(self.tabix_cmd)


class Binning(SGEJobTask):
    """
    Move to annodb pipeline 
    Task to run the coverage binning script
    """


    bam_file = luigi.Parameter(default='bam')
    gvcf = luigi.Parameter(default='gvcf')
    sample = luigi.Parameter()
    block_size = luigi.Parameter()
    scratch = luigi.Parameter()
    mode = luigi.Parameter()
    
    
    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/rp2801/genome_cvg_temp"
    
    
    def __init__(self,*args,**kwargs):

        super(Binning,self).__init__(*args,**kwargs)

        if not os.path.isdir(self.scratch): ## Recursively create the directory if it doesnt exist
            os.makedirs(self.scratch)

        self.genomecov_bed = os.path.join(self.scratch,self.sample+'.genomecvg.bed')
        self.human_chromosomes = []
        self.human_chromosomes.extend(range(1, 23))
        self.human_chromosomes = [str(x) for x in self.human_chromosomes]
        self.human_chromosomes.extend(['X', 'Y','MT'])

        if self.mode == 'coverage':
            self.binning_cmd = "{0} {1} --sample_id {2} --block_size {3} --output_dir {4} {5} --mode coverage".format(
                config().pypy_loc,config().binner_loc,self.sample,self.block_size,self.scratch,
                self.genomecov_bed)        
        elif self.mode == 'gq':
            self.binning_cmd = "{0} {1} --sample_id {2} --block_size {3} --output_dir {4} {5} --mode gq".format(
                config().pypy_loc,config().binner_loc,self.sample,self.block_size,self.scratch,
                self.gvcf)
   
    def requires(self):
        """
        Dependency is the completion of the CreateGenomeBed Task
        """
    
        if self.mode == 'coverage':
            return [CreateGenomeBed(bam=self.bam_file,sample=self.sample,scratch=self.scratch)]
        elif self.mode == 'gq':
            return MyExtTask(self.gvcf)
        
    
    def output(self):
        """
        Output from this task are 23 files for each chromosome
        """

        
        for chrom in self.human_chromosomes:
            file_loc = os.path.join(self.scratch,self.sample+'_%s_binned_'%self.mode+self.block_size
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
    

    bam = luigi.Parameter()
    sample_name = luigi.Parameter()
    scratch = luigi.Parameter()
    
    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"

    
    def __init__(self,*args,**kwargs):
        super(AlignmentMetrics,self).__init__(*args,**kwargs)

        #self.sample_dir = os.path.join(self.scratch,self.sample_name)
        

        if not os.path.isdir(self.scratch): ## Recursively create the directory if it doesnt exist
            os.makedirs(self.scratch)

        self.output_file_raw = os.path.join(self.scratch,self.sample_name+'.alignment.metrics.raw.txt')
        self.output_file = os.path.join(self.scratch,self.sample_name+'.alignment.metrics.txt')
        self.cmd = "{0} -XX:ParallelGCThreads=12 -jar {1} CollectAlignmentSummaryMetrics TMP_DIR={2} VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE={3} INPUT={4} OUTPUT={5}".format(config().java,config().picard,self.scratch,config().ref,self.bam,self.output_file_raw)               
        self.parser_cmd = """cat {0} | grep -v "^#" | awk -f /home/rp2801/git/dragen_pipe/dragen/python/transpose.awk > {1}""".format(self.output_file_raw,self.output_file)
        
        
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
        os.system(self.parser_cmd)
        

class DuplicateMetrics(SGEJobTask):
    """
    Parse Duplicate Metrics from dragen logs
    """
    
    dragen_log = luigi.Parameter()
    sample_name = luigi.Parameter()
    

    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch11/rp2801/genome_cvg_temp"

    
    def __init__(self):
        super(DuplicateMetrics,self).__init__(*args,**kwargs)
        if not os.path.isdir(self.scratch): ## Recursively create the directory if it doesnt exist
            os.makedirs(self.scratch)

        self.output_file = os.path.join(self.scratch,self.sample_name+'duplicates.txt')
        self.cmd = "grep 'duplicates marked' {0}".format(self.dragen_log)
        

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
        if match:
            perc_duplicates = match.group(1)
            with open(self.output_file,'w') as OUT_HANDLE:
                print >> OUT_HANDLE,perc_duplicates
                
        else: ## Implement logging
            print "Could not find duplicate metrics in dragen log!"

        
class VariantCallingMetrics(SGEJobTask):
    """
    Run the picard tool for qc evaluation of the variant calls 
    """

    sample_name = luigi.Parameter()
    vcf = luigi.Parameter()
    scratch = luigi.Parameter()

    ## System Parameters 
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"
      
    
    def __init__(self,*args,**kwargs):
        super(VariantCallingMetrics,self).__init__(*args,**kwargs)
        #sample_name = luigi.Parameter()
        if not os.path.isdir(self.scratch): ## Recursively create the directory if it doesnt exist
            os.makedirs(self.scratch)

        self.output_file_raw = os.path.join(self.scratch,"{0}.variant_calling_summary_metrics.raw.txt".format(self.sample_name))
        self.output_file = os.path.join(self.scratch,"{0}.variant_calling_summary_metrics.txt".format(self.sample_name))
        self.cmd = "{0} -XX:ParallelGCThreads=12 -jar {1} CollectVariantCallingMetrics INPUT={2} OUTPUT={3} DBSNP={4}".format(config().java,config().picard,self.vcf,self.output_file_raw,config().dbsnp)
        self.parser_cmd = """cat {0} | grep -v '^#'| awk -f /home/rp2801/git/dragen_pipe/dragen/python/transpose.awk > {1} """.format(self.output_file_raw,self.output_file)

        
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
        
        return luigi.LocalTarget("{0}.processed".format(self.output_file))

    def work(self):
        """
        Run this task
        """
    
        os.system(self.cmd)
        os.system(self.parser_cmd) 


class CheckUpdates(luigi.Target):
    """
    Checks whether the database entry was updated
    correctly
    """

    
    def __init__(self,*args,**kwargs):
        super(CheckUpdates,self).__init__(*args,**kwargs)
        cnf_file = kwargs['cnf_file']
        value = kwargs['value']
        db_field = kwargs['db_field']
        seq_type = kwargs['seq_type']
        prep_id = kwargs['prep_id']
        chgvid = kwargs['chgvid']
        testmode = kwargs['testmode']
        
        
        if self.testmode:
            db = MySQLdb.connect(read_default_group='clienttestdb',db='sequenceDB',read_default_file=self.cnf_file)
            db_cursor = db.cursor()
            
        query_statement = """ SELECT {0} FROM seqdbClone WHERE CHGVID = '{1}' AND seq_type = '{2}' AND prepid = {3} """.format(self.db_field,self.chgvid,self.seq_type,self.prep_id)
        db_cursor.execute(query_statement)
        self.db_val = testdb_cursor.fetchall()
        
    
    def exists(self):
        if len(db_val > 1):
            print "More than 1 entry , warning : duplicate prepids !"
        if self.value == self.db_val[0][-1]:
            return True
        else:
            return False
        
    
class UpdateDatabase(SGEJobTask):
    """
    Populate database with output files
    """

    output_file = luigi.Parameter()
    ## A dict parameter containing the field names to be used for update
    ## from the output file and the corresponding field in the database 
    fields = luigi.DictParameter() 
    cnf_file = luigi.Parameter()
    testmode = luigi.BoolParameter()
    seqtype = luigi.Parameter()
    prepid = luigi.Parameter()
    chgvid = luigi.Parameter()  
        
    
    def __init__(self,*args,**kwargs):
        super(UpdateDatabase,self).__init__(*args,**kwargs)
        if self.testmode:
            db = MySQLdb.connect(read_default_group='clienttestdb',db='sequenceDB',read_default_file=self.cnf_file)
            db_cursor = db.cursor()

        self.query_statement = """ UPDATE seqdbClone SET {0} = {1} WHERE CHGVID = {2} AND seqType = {3} AND prepid = {4}"""
        
        self.out_hash = {}
        with open(self.output_file,'w') as OUT:
            for line in OUT:
                line = line.strip('\n')
                contents = line.split(' ')
                self.out_hash[contents[0]] = contents[1]


    def requires(self):
        """
        The requirement for this task
        the outputfile from different tasks should
        be present
        """
        
        return [MyExtTask(self.output_file)]
    

    def output(self):
        """
        The output from this task
        """

        check_args = {'cnf_file':self.cnf_file,'seq_type':self.seq_type,
                      'prep_id':self.prep_id,'chgvid':self.chgvid,
                      'testmode':self.testmode}
        for key in self.fields.keys():
            check_args['db_field'] = self.fields[key]
            check_args['value'] = self.out_hash[key]
            return CheckUpdates(check_args)
 
    def work(self):
        """
        Run this task
        """


        for key in self.fields.keys():
            value = self.out_hash[field]
            db_field = self.fields[key]
            self.db_cursor.execute(self.query_statement.format(db_field,value,self.chgvid,self.seqtype,self.prepid))
            self.db.commit()
            
        self.db.close()

       
class GenderCheck(SGEJobTask):
    """
    Check gender using X,Y chromosome coverage and update seqgender field
    """

    ## Parameters for this task
    sample = luigi.Parameter()
    seqtype = luigi.Parameter()
    prepid = luigi.Parameter()
    testmode = luigi.BoolParameter()
    ## Contains defaults for dragen pipeline
    
    
    def __init__(self,*args,**kwargs):
        super(GenderCheck,self).__init__(*args,**kwargs)
        ## Possible set of update entries for the seqgender column
        self.valid_updates = ['M','F','Ambiguous','N/A']
        
        ## Read the cnf_file and get a connection to the database
        if self.testmode:
            self.db = MySQLdb.connect(read_default_group='clienttestdb',db='sequenceDB',read_default_file=self.cnf_file)
            self.db_cursor = db.cursor()
            
        ## Define the SQL queries that need to be executed in this task
        
        ## Get X,Y Coverages
        temp_query = ("""SELECT XCvg, YCvg FROM seqdbClone WHERE CHGVID = {0}
                      AND seqType = {1} AND prepID = {2}"""
                     )
        
        self.db_cursor.execute(temp_query)
        x_cvg = query_result[0][1]
        y_cvg = query_result[0][2]
        
        sample_query = ("""SELECT CHGVID FROM seqdbClone WHERE CHGVID = {0} AND
                        AND seqType = {1} AND prepID = {2}"""
                        )
    
    
        
    def requires(self):
        """
        The requirement for this task 
        """
        ## Define a better dependency
        #self.

        ## Run Cvg Metrics on X and Y chromosomes
        
        
        
    def output(self):
        """
        The output from this task
        """

        ## The sucess of the task depends on the sample seqgender being
        ## updated to one of the valid set of update entries (possible
        ## potential for a bug here where a wrong update could be made and
        ## go unoticed in this check)

        ## Get a cursor to the database
        temp_cursor = self.db.cursor()
        
        
    def work(self):
        """
        Run this task
        """

        
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
    ped_loc = luigi.Parameter()

    
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"

    
    def __init__(self,*args,**kwargs):
        super(AppendMasterPed,self).__init__(*args,**kwargs)
        self.cmd = "cat {0} >> {1}"
        
    def requires(self):
        """
        Make sure relevant files are present
        """

        return [MyExtTask(self.masterped),MyExtTask(self.ped_loc)]

    def work(self):
        """
        Run the task
        """
        os.system(self.cmd.format(self.ped_loc,self.masterped))
        
    
class RunRelatednessCheck(SGEJobTask):
    """
    This Task would be run nightly as a batch 
    """

    masterped = luigi.Parameter()
    threshold_file = luigi.Parameter()
    relatedness_map = luigi.Parameter()
    database = luigi.Parameter()
    cyrptic_threshold = luigi.Parameter()
    cnf = luigi.Parameter() ## MySQL config file 
    database_table = luigi.Parameter()
    test_mode = luigi.Parameter()
    #last_run = luigi.DateParameter()
    
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"

    
    def __init__(self,*args,**kwargs):
        super(RunRelatednessCheck,self).__init__(*args,**kwargs)
        self.append_cmd = "cat {0} >> {1}"
        self.relatedness_cmd = "{0} new_relatedness_check.py {1} {2} {3} {4} {5}".format(config().python_path,config().scratch,self.threshold_file,self.cryptic_threshold,self.test_mode)
        ## Check the database for the number of samples to be appended
        db = MySQLdb.connect(db=self.database,read_default_group="clientseqdb",read_default_file=self.cnf)
        ## Get the samples to append to master ped 
        db_cursor = db.cursor()        
        temp_query = "USE %s SELECT CHGVID,ped_location FROM %s WHERE is_append = False"%(database,database_table)
        db_cursor.execute(temp_query)
        self.result = db_cursor.fetchall()
                

        
    def requires(self):
        """
        The dependencies for this task
        """

        for chgvid,ped_loc in self.result:
            yield AppendMasterPed(masterped,ped_loc)
            
        yield MyExtTask(threshold_file)
        yield MyExtTask(relatedness_map)
        
    def work(self):
        """
        Run the python relatedness script 
        """
        
        ## Iterate over the results and append the individual peds to the masterped
        for element in self.result:
            ped_location = element[1]
            os.system(self.append_cmd.format(ped_location,self.masterped))

        ## Run the relatedness pipeline
        os.system(self.relatedness_cmd)
            
    def output(self):
        """
        The output from this task
        The line count should be appended by 1 
        """

        return 


class CheckAppend(luigi.Target):
    """
    A Target class for checking if the append was sucessful
    """
    
    masterped = luigi.Parameter()
    num_samples = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"

    def __init__(self,*args,**kwargs):
        super(CheckAppend,self).__init__(*args,**kwargs)
        
