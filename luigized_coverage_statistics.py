#!/usr/bin/env Python
import sys
import subprocess
import os
import utils 
import luigi
from luigi.util import inherits 
from luigi.contrib.sge import SGEJobTask,LocalSGEJobTask
import logging
import time

##################################################################################################
##        Run as :                                                                              ##  
##        luigi --module luigized_coverage_statistics RunCvgMetrics  --local-scheduler          ##
##        Parameters are stored in luigi.cfg                                                    ##
##                                                                                              ##
##################################################################################################

## Todo :
## 1. Add new tasks from redmine
## 2. Change os to subprocess
## 3. Look into how root task is initialized in luigi
## 4. Define better dependencies 
## 5. To integrate with GATK pipeline , look into http://datapipelinearchitect.com/luigi-non-input-task-dependencies/
##    The idea is to have tasks which need to be completed without worrying about checking for output dependencies

################################################################################################################################################


#logger = logging.getLogger('luigi-interface')



class config(luigi.Config):
    """
    config class for instantiating parameters for this pipeline
    the values are read from luigi.cfg in the current folder 
    """

    java = luigi.Parameter(default='/nfs/goldstein/software/jdk1.8.0_05/bin/java')
    picard = luigi.Parameter(default='/nfs/goldstein/software/picard-tools-2.2.1/picard.jar')
    reference_file = luigi.Parameter(default='/nfs/goldsteindata/refDB/HS_Build37/BWA_INDEX_hs37d5/hs37d5.fa')
    seqdict_file = luigi.Parameter(default='/nfs/goldsteindata/refDB/HS_Build37/BWA_INDEX_hs37d5/hs37d5.dict')
    bed_file = luigi.Parameter(default='=/home/rp2801/test_folder/coverage_statistics/ccds_r14.bed')
    target_file = luigi.Parameter(default='/home/rp2801/test_folder/coverage_statistics/CCDS.genes.unique.regions.list')
    create_targetfile = luigi.BooleanParameter(default=False)
    bait_file = luigi.Parameter(default='/home/rp2801/test_folder/coverage_statistics/SeqCap_EZ_Exome_v3_capture.bed.list')
    wgsinterval = luigi.BooleanParameter(default=False)
    scratch = luigi.Parameter(default='./')
    output_dir = luigi.Parameter(default='./')
    sample = luigi.Parameter(default='Schiz6113363')
    bam_file = luigi.Parameter(default='/nfs/seqscratch11/filtered.Schiz6113363.production.bam')
    seq_type = luigi.Parameter(default='Genome')
    python_path = luigi.Parameter()
    relatedness_refs = luigi.Parameter()
    sampleped_loc = luigi.Parameter()
    relatedness_markers = luigi.Parameter()
    bedtools_loc = luigi.Parameter()
    pypy_loc = luigi.Parameter()
    coverage_binner_loc = luigi.Parameter()
    

class RootTask(luigi.WrapperTask):
    """
    Wrapper Task
    """

    ## Read in the bam file from the cmd line just to make sure
    ## to ensure uniqueness of the RootTask 
    bam_file = luigi.Parameter()
    
    
    def requires(self):
        self.output_file = config().scratch+config().sample+'.cvg.metrics.txt'
        return [ParsePicardOutputs(self.output_file,self.bam_file)]



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


    sample = luigi.Parameter()
    bam_file = luigi.Parameter()
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
        self.cvg_file = "{scratch}{samp}.cvg.metrics.txt".format(
            scratch=config().scratch, samp=self.sample)

        self.output_file = "{scratch}{samp}.cvg.metrics.txt.transposed".format(
            scratch=config().scratch, samp=self.sample)

        ## Define the command to be run 
        self.cmd = """cat %s | grep -v '^#'| head -3 | awk -f /home/rp2801/git/dragen_pipe/transpose.awk > %s """%(self.output_file,self.output_file)

    def requires(self):
        """ 
        Dependency is just the presence of the awk script
        and the completion of the RunCvgMetrics Task
        """
        
        return RunCvgMetrics(self.bam_file,self.sample,self.seq_type)

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
    
    
    bam_file = luigi.Parameter()
    sample = luigi.Parameter()
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
        self.output_file = "{scratch}{samp}.cvg.metrics.txt".format(
            scratch=config().scratch, samp=self.sample)

    
        ## Initialize the targetfile name if it needs to be created
        if config().create_targetfile == True:
            self.target_file = ((config().scratch) +
                           (utils.get_filename_from_fullpath(config().bed_file)) +
                           (".list"))
        else:
            self.target_file = config().target_file 

        ## Setup the command to be run 

        if self.seq_type.upper() == 'GENOME': ## If it is a genome sample
            if self.wgsinterval == True: ## Restrict wgs metrics to an interval
                cvg_cmd = ("%s -jar %s CollectWgsMetrics R=%s O=%s I=%s INTERVALS=%s MQ=20 Q=10"%(config().java,config().picard,config().reference_file,self.output_file,self.bam_file,self.target_file))
                
            else: ## Run wgs metrics across the genome
                cvg_cmd = ("%s -jar %s CollectWgsMetrics R=%s O=%s I=%s MQ=20 Q=10"%(config().java,config().picard,config().reference_file,self.output_file,self.bam_file))
              
                             
        else: ## Anything other than genome i.e. Exome or Custom Capture( maybe have an elif here to check further ?)
            cvg_cmd = ("%s -jar %s CalculateHsMetrics BI=%s TI=%s METRIC_ACCUMULATION_LEVEL=ALL_READS I=%s O=%s MQ=20 Q=10"%(config().java,config().picard,config().bait_file,self.target_file,self.bam_file,self.output_file))

        
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
            yield [CreateTargetFile(),MyExtTask(self.bam_file)]
        else:
            yield MyExtTask(self.bam_file)
        
        
    def work(self):
        """
        Run Picard CalculateHsMetrics or WgsMetrics
        """
        ## Run the command
        os.system(cvg_cmd)
        #os.system("touch %s"%(self.output_file))


class CreateTargetFile(SGEJobTask):
    """ 
    A Luigi Task
    """
    

    def __init__(self,*args,**kwargs):
        super(CreateTargetFile,self).__init__(*args,**kwargs)

        temp = utils.get_filename_from_fullpath(config().bed_file) + ".list"
        self.output_file = "{scratch}{temp_name}.genomecvg.bed".format(
                            scratch=config().scratch, temp_name=temp)

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
        os.system(bedtointerval_cmd)

        
class CreatePed(SGEJobTask):
    """ 
    A luigi Task for creating ped files
    """
    
    ped_type = luigi.Parameter() ## This has one of two values :
                                 ## ethnicity or relatedness
    bam_file = luigi.Parameter()
    vcf_file = luigi.Parameter()
    sample = luigi.Parameter()
    family_id = luigi.Parameter(default="")
    seq_type = luigi.Parameter()
    
    
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"


    def __init__(self,*args,**kwargs):
        super(CreatePed,self).__init__(*args,**kwargs)

        self.output_file = config().scratch+self.sample+'.'+self.ped_type+'.ped'

        if self.ped_type == 'relatedness':
            self.sampletoped_cmd = ("%s %s --scratch %s --vcffile %s --bamfile %s" 
                               " --markers %s --refs %s --sample %s --fid %s"
                               " --seq_type %s --ped_type %s --output %s --pythondir /home/rp2801/git/dragen_pipe/"
                               %(config().python_path,config().sampleped_loc,
                                 config().scratch,self.vcf_file,
                                 self.bam_file,config().relatedness_markers,
                                 config().relatedness_refs,self.sample,
                                 self.family_id,self.seq_type,self.ped_type,
                                 config().scratch)
                               )
            
        elif self.ped_type == 'ethnicity':
            self.sampletoped_cmd = ("%s %s --scratch %s --vcffile %s --bamfile %s"
                               " --markers %s --refs %s --sample %s"
                               " --seq_type %s --ped_type %s --output %s" 
                               %(config().python_path,config().sampleped_loc,
                               config().scratch,self.vcf_file,
                               self.bam_file,config().relatedness_markers,
                               config().relatedness_refs,self.sample,
                               self.seq_type,self.ped_type,config().scratch)
                               )


    def requires(self):
        """
        The dependencies for this task are the existence of the relevant files
        """
        
        return [MyExtTask(self.bam_file),MyExtTask(self.vcf_file),MyExtTask(config().relatedness_refs)]


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

    bam_file = luigi.Parameter()
    sample = luigi.Parameter()
    
    def __init__(self,*args,**kwargs):
        super(CreateGenomeBed,self).__init__(*args,**kwargs)

        self.genomecov_bed = "{scratch}{samp}.genomecvg.bed".format(
                              scratch=config().scratch, samp=self.sample)
        self.genomecov_cmd = "%s genomecov -bga -ibam %s > %s"%(
                              config().bedtools_loc,self.bam_file,
                              self.genomecov_bed)    
   
      
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"

    
    def requires(self):
        """
        The dependency is the presence of the bam file
        """
        
        return MyExtTask(self.bam_file)

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

        self.genomecov_bed = "{scratch}{samp}.genomecvg.bed".format(
            scratch=config().scratch, samp=self.sample)

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
       
        return [CreateGenomeBed(self.bam_file,self.sample)]
        
    
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

        print self.binning_cmd
        os.system(self.binning_cmd)


class AlignmentMetrics(LocalSGEJobTask):
    """
    Parse Alignment Metrics from dragen logs
    """
    

    bam_file = luigi.Parameter()
    sample = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"

    def __init__(self):
        super(AlignmentMetrics,self).__init__(*args,**kwargs)
        self.output_file = '{scratch}{samp}.alignment.metrics.txt'.format(scratch=config().scratch,samp=self.sample)
        self.cmd = "%s %s CollectAlignmentSummaryMetrics TMP_DIR=%s VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=%s INPUT=%s OUTPUT=%s"%(config().java,config().picard,config().scratch,config().reference_file,self.bam_file,self.output_file)               
    
    
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
        
        
        return MyExtTask(self.bam_file) 

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
    sample = luigi.Parameter()
    

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"

    
    def __init__(self):
        super(AlignmentMetrics,self).__init__(*args,**kwargs)
        #self.output_file = '{scratch}{samp}.alignment.metrics.txt'.format(scratch=config().scratch,samp=self.sample)
        self.cmd = "grep 'duplicates marked' %s"%self.dragen_log
        
    def requires(self):
        """ 
        Dependencies for this task is the existence of the log file 
        """

        return MyExtTask(self.bam_file)


    def output(self):
        """
        """

        return luigi.LocalTarget('dump.txt')

    def work(self):
        """
        Execute this task
        """
        ## The regular expression to catch the percentage duplicates in the grep string
        catch_dup = re.compile('.*\((.*)%\).*')

        dragen_output = subprocess.check_output(self.cmd,shell=True)
        match = re.match(catch_dup,dragen_output)
        perc_duplicates = match.group(1)
        

def EthnicityPipeline(LocalSGEJobTask):
    """
    Run the Ethnicity Pipeline
    """


class CheckAppend(luigi.Target):
    """
    This is a target class for checking whether the database
    entry was correctly flagged
    """


    sample = luigi.Parameter()
    seq_type = luigi.Parameter()


    def __init__(self):
        super(AppendMasterPed,self).__init__(*args,**kwargs)
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
class AppendMasterPed(SGEJobTask):
    """
    Append the relatedness ped to the master ped,
    A DB entry is created stating that the sample needs
    to have the ped file appended 
    The actual append is performed as part of the relatedness
    pipeline 
    """

    ped_type = luigi.Parameter(significant=True) ## This has one of two values :
                                 ## ethnicity or relatedness
    bam_file = luigi.Parameter(significant=True)
    vcf_file = luigi.Parameter(significant=True)
    sample = luigi.Parameter(significant=True)
    family_id = luigi.Parameter(default="",significant=True)
    seq_type = luigi.Parameter(significant=True)
    master_ped = luigi.Parameter()


    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"
      
    
    def __init__(self,*args,**kwargs):
        super(AppendMasterPed,self).__init__(*args,**kwargs)       
        self.output_file = config().scratch+self.sample+'.'+self.ped_type+'.ped'        
        self.cmd = '''seqdb -e "use sequenceDB; UPDATE table_name SET to_append=TRUE WHERE sample_name=%s, seq_type=%s" '''%(self.sample,self.seq_type)


    def requires(self):
        """
        The dependency for this task is the completion of the
        CreatePed Task with the appropriate signatures
        """
        
        return CreatePed(self.ped_type,self.bam_file,self.vcf_file,
                         self.sample,self.family_id,self.seq_type)

    
    def output(self):
        """
        The Output from this task 
        Define the target,such that
        the line count has increased by 1
        in the masterped file 
        """
        
        return CheckAppend(self.sample,self.seq_type)
        

    def work(self):
        """
        Run this task, a simple append
        should do it
        """

        os.system(self.cmd)

        

class VariantCallingMetrics(SGEJobTask):
    """
    Run the picard tool for qc evaluation of the variant calls 
    """

    vcf_file = luigi.Parameter()
    sample = luigi.Parameter()

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"
      
    
    def __init__(self,*args,**kwargs):
        super(AppendMasterPed,self).__init__(*args,**kwargs)
        self.output_file = self.sample + '.variant.metrics.txt'
        self.cmd = "java -jar %s CollectVariantCallingMetrics INPUT=%s OUTPUT=%s DBSNP=%s"%(config().picard,self.vcf_file,config().dbsnp)

    
    def requires(self):
        """
        The requirement for this task is the presence of the analysis ready 
        vcf file from the GATK pipeline
        """
    
        return MyExtTask(self.vcf_file)

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
    

class UpdateDatabase(SGEJobTask):
    """
    Populate database with output files
    """

    self.output_file = luigi.Parameter()

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

        
