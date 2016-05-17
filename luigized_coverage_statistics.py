#!/usr/bin/env Python
import sys
import subprocess
import os
import utils 
import luigi

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



class config(luigi.Config):
    java = luigi.Parameter()
    picard = luigi.Parameter()
    reference_file = luigi.Parameter()
    seqdict_file = luigi.Parameter()
    bed_file = luigi.Parameter()
    bam_file = luigi.Parameter()
    target_file = luigi.Parameter()
    create_targetfile = luigi.BooleanParameter()
    seq_type = luigi.Parameter()
    bait_file = luigi.Parameter()
    sample = luigi.Parameter()
    wgsinterval = luigi.BooleanParameter()
    scratch = luigi.Parameter()
    output_dir = luigi.Parameter()


class MyExtTask(luigi.ExternalTask):

    file_loc = luigi.Parameter()
    def output(self):
        return luigi.LocalTarget(self.file_loc)

    
class RunCvgMetrics(luigi.Task):
    """ 
    A luigi task
    """
    
    
    output_file = config().scratch+config().sample+'cvg.metrics.txt'
    if config().create_targetfile == True:
        target_file = ((config().scratch) +
        (utils.get_filename_from_fullpath(config().bed_file)) +
        (".list"))
    else:
        target_file = config().target_file
        
        
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

        if config().wgsinterval == True:
            return MyExtTask(config().reference_file)
        else:
            return CreateTargetFile()
                
        
    def run(self):
        """
        Run Picard CalculateHsMetrics or WgsMetrics
        """
        
        if config().seq_type.upper() == 'GENOME': ## If it is a genome sample
            if config().wgsinterval == True: ## Restrict wgs metrics to an interval
                wgs_cmd = ("%s -jar %s CollectWgsMetrics R=%s O=%s I=%s"
                           "INTERVALS=%s MQ=20 Q=10"
                           %(config().java,
                             config(),picard,
                             config().reference_file,
                             config().output_file,
                             config().bam_file,
                             self.target_file))
                
            else: ## Run wgs metrics across the genome
                wgs_cmd = ("%s -jar %s CollectWgsMetrics R=%s O=%s I=%s"
                           "MQ=20 Q=10"
                           %(config().java,
                             config().picard,
                             config().reference_file,
                             config().output_file,
                             config().bam_file))

                
            #os.system(wgsmetrics_cmd)
            os.system("touch %s"%(self.output_file))
                             
        else: ## Anything other than genome i.e. Exome or Custom Capture 
            hsmetrics_cmd = ("%s -jar %s CalculateHsMetrics BI=%s TI=%s" 
                         "METRIC_ACCUMULATION_LEVEL=ALL_READS I=%s O=%s"
                         "MQ=20 Q=10"%(config().java,
                                       config().picard,
                                       config().bait_file,
                                       self.target_file,
                                       config().bam_file,
                                       self.output_file))
            ## Run the command
            #os.system(hsmetrics_cmd)
            os.system("touch %s"%(self.output_file))


class CreateTargetFile(luigi.Task):
    """ 
    A Luigi Task
    """
      
    output_file = (config().scratch +
                   utils.get_filename_from_fullpath(config().bed_file) +
                   ".list")
    
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
       
        bedtointerval_cmd = ("%s -jar %s BedToIntervalList I=%s SD=%s" 
                             "OUTPUT=%s"%(config().java,
                                         config().picard,
                                         config().bed_file,
                                         config().seqdict_file,
                                         self.output_file))    
        
        os.system("touch %s"%self.output_file)
        #os.system(bedtointerval_cmd)

        
    

