#!/usr/bin/env Python
import sys
import subprocess
import os
import utils 
import luigi
from luigi.util import inherits 
from luigi.contrib.sge import SGEJobTask,LocalSGEJobTask
import logging

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


class ParsePicardOutputs(LocalSGEJobTask):
    """
    Simple Parser for Picard output
    """


    output_file = luigi.Parameter()
    bam_file = luigi.Parameter()
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"


    def requires(self):
        """ 
        Dependency is just the presence of the awk script
        and the completion of the RunCvgMetrics Task
        """
        
        return RunCvgMetrics(self.bam_file,self.output_file)

    def output(self):
        """
        Output from this task
        """

        return luigi.LocalTarget(self.output_file+'.transposed')

    def work(self):
        """
        Just a simple shell command is enough!
        """
        
        cmd = """cat %s | grep -v '^#'| head -3 | awk -f /home/rp2801/git/dragen_pipe/transpose.awk > %s.transposed """%(self.output_file,self.output_file)
        #os.system(cmd)
        os.system("""touch %s.transposed"""%(self.output_file))


class RunCvgMetrics(LocalSGEJobTask):
    """ 
    A luigi task
    """
    
    
    bam_file = luigi.Parameter(default='/nfs/seqscratch11/filtered.Schiz6113363.production.bam')
    #output_file = config().scratch+config().sample+'.cvg.metrics.txt'
    output_file = luigi.Parameter()
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/home/rp2801/git"

        
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
            self.target_file = ((config().scratch) +
                           (utils.get_filename_from_fullpath(config().bed_file)) +
                           (".list"))
        else:
            self.target_file = config().target_file


        if config().create_targetfile == True:
            yield [CreateTargetFile(),MyExtTask(self.bam_file)]
        else:
            yield MyExtTask(self.bam_file)
        
        
    def work(self):
        """
        Run Picard CalculateHsMetrics or WgsMetrics
        """
 
        print config().seq_type
        print config().create_targetfile 

        if config().seq_type.upper() == 'GENOME': ## If it is a genome sample
            if config().wgsinterval == True: ## Restrict wgs metrics to an interval
                wgs_cmd = ("%s -jar %s CollectWgsMetrics R=%s O=%s I=%s INTERVALS=%s MQ=20 Q=10"%(config().java,config().picard,config().reference_file,self.output_file,self.bam_file,self.target_file))
                
            else: ## Run wgs metrics across the genome
                wgs_cmd = ("%s -jar %s CollectWgsMetrics R=%s O=%s I=%s MQ=20 Q=10"%(config().java,config().picard,config().reference_file,self.output_file,self.bam_file))

                
            #os.system(wgs_cmd)
            os.system("touch %s"%(self.output_file))
                             
        else: ## Anything other than genome i.e. Exome or Custom Capture( maybe have an elif here to check further ?)
            hsmetrics_cmd = ("%s -jar %s CalculateHsMetrics BI=%s TI=%s METRIC_ACCUMULATION_LEVEL=ALL_READS I=%s O=%s MQ=20 Q=10"%(config().java,config().picard,config().bait_file,self.target_file,self.bam_file,self.output_file))
            ## Run the command
            os.system(hsmetrics_cmd)
            #os.system("touch %s"%(self.output_file))


class CreateTargetFile(SGEJobTask):
    """ 
    A Luigi Task
    """
      


    output_file = (config().scratch +
                   utils.get_filename_from_fullpath(config().bed_file) +
                   ".list")

    #output_file = 'chutiye!!!'
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

    def run(self):
        """
        Run Picard BedToIntervalList 
        """
       
        bedtointerval_cmd = ("%s -jar %s BedToIntervalList I=%s SD=%s OUTPUT=%s"%(config().java,config().picard,config().bed_file,config().seqdict_file,self.output_file))    
        
        os.system("touch %s"%self.output_file)
        #os.system(bedtointerval_cmd)

        
class CreatePed(LocalSGEJobTask):
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
    output_dir = luigi.Parameter(default="./")
    

    def requires(self):
        """
        The dependencies for this task are the existence of the relevant files
        """
        
        return [MyExtTask(self.bam_file),MyExtTask(self.vcf_file),MyExtTask(config().relatedness_refs)]

    def output(self):
        """
        Output from this task is a ped file named with correct type
        """
        
    
        return luigi.LocalTarget(config().scratch+self.sample+'.'+self.ped_type+'.ped')

    def work(self):
        """
        To run this task, the sampletoped_fixed.py script is called
        """

        if self.ped_type == 'relatedness':
            sampletoped_cmd = ("%s %s --scratch %s --vcffile %s --bamfile %s" 
                               " --markers %s --refs %s --sample %s --fid %s"
                               " --seq_type %s --ped_type %s --output %s"
                               %(config().python_path,config().sampleped_loc,
                                 config().scratch,self.vcf_file,
                                 self.bam_file,config().relatedness_markers,
                                 config().relatedness_refs,self.sample,
                                 self.family_id,self.seq_type,self.ped_type,
                                 config().scratch)
                               )
            
        elif self.ped_type == 'ethnicity':
            sampletoped_cmd = ("%s %s --scratch %s --vcffile %s --bamfile %s"
                               " --markers %s --refs %s --sample %s"
                               " --seq_type %s --ped_type %s --output %s" 
                               %(config().python_path,config().sampleped_loc,
                               config().scratch,self.vcf_file,
                               self.bam_file,config().relatedness_markers,
                               config().relatedness_refs,self.sample,
                               self.seq_type,self.ped_type,config().scratch)
                               )
            
        os.system(sampletoped_cmd)
        #os.system("%s sampletoped_fixed.py --scratch %s")

        
class CreateGenomeBed(LocalSGEJobTask):
    """
    Use bedtools to create the input bed file for the coverage binning script
    """

    bam_file = luigi.Parameter()
    sample = luigi.Parameter()
    genomecov_bed = luigi.Parameter()
    
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

        genomecov_cmd = "%s genomcov -bga -ibam %s > %s"%(config().bedtools_loc,self.bam_file,
                                                          self.genomecov_bed)

        os.system(genomecov_cmd)


class CoverageBinning(LocalSGEJobTask):
    """
    Task to run the coverage binning script
    """
    bam_file = luigi.Parameter()
    sample = luigi.Parameter()
    block_size = luigi.Parameter()
    

    def __init__(self,*args,**kwargs):
        self.genomecov_bed = config().scratch()+self.sample+'.genome.cvg'
        self.human_chromosomes.extend(range(1, 23))
        self.human_chromosomes = [str(x) for x in human_chromosomes]
        self.human_chromosomes.extend(['X', 'Y','MT'])
        
    def requires(self):
        """
        Dependency is the completion of the CreateGenomeBed Task
        """
       
        return CreateGenomeBed(self.bam_file,self.sample,self.genomecov_bed)

    
    def output(self):
        """
        Output from this task are 23 files for each chromosome
        """

        for chrom in self.human_chromosomes:
            yield [luigi.LocalTarget(self.sample+"_read_coverage_" +
                                     self.block_size + "_chr%s.txt"%chrom)]

    def work(self):
        """ Run the binning script
        """

        binning_cmd = "%s %s %s %s"%(config().pypy_loc,
                                     config().coverage_binner_loc,
                                     self.block_size,self.genomecov_bed)
