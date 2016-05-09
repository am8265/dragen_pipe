#!/usr/bin/env Python
import sys
import subprocess
import os
import utils 
import luigi

##################################################################################################
##        Run as :                                                                              ##  
##        luigi --module luigized_coverage_statistics RootTask --args                           ##
##        '{"java":"java location", "picard":"picard location", "bed_file":"bed location"...}'  ##
##        --local-scheduler                                                                     ##
##################################################################################################

## Todo :
## 1. Add new tasks from redmine
## 2. Change os to subprocess
## 3. Look into how root task is initialized in luigi





class CreateTargetFile(luigi.Task):
    """ 
    A Luigi Task
    """

    args = luigi.DictParameter()            
    output_file = (self.args['output_dir'] +
                   utils.get_filename_from_fullpath(self.args['bed_file']) +
                   ".list")
    
    def requires(self):
        """ 
        Dependcy for this task is the existence of the bedfile
        and the human build37 sequence dictionary file
        """
        
        return [luigi.LocalTarget(self.args['bed_file']),luigi.LocalTarget(self.args['seqdict_file'])]
    
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
                             "OUTPUT=%s"%(self.args['java_location'],
                                         self.args['picard_location'],
                                         self.args['bed_file'],
                                         self.args['seqdict_file'],
                                         self.output_file))    
        
        os.system("touch %s"%self.output_file)
        #os.system(bedtointerval_cmd)

        
class RunMetrics(luigi.Task):
    """ 
    A luigi task
    """

    args = luigi.DictParameter()
    output_file = args['scratch'] + args['sample'] + '.hsmetrics.txt'
    
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
        if self.args['create_targetfile'] == True:
            return CreateTargetFile()
        
        if self.args['wgsinterval'] == True:
            return luigi.LocalTarget(self.args['reference_file'])
    
    def run(self):
        """
        Run Picard CalculateHsMetrics or WgsMetrics
        """
        
        if self.args['seq_type'].upper() == 'GENOME':
            if self.args['wgsinterval'] == True: ##
                wgs_cmd = ("%s -jar %s CollectWgsMetrics R=%s O=%s I=%s"
                           "INTERVALS=%s MQ=20 Q=10"
                           %(self.args['java_location'],
                             self.args['picard_location'],
                             self.args['reference_file'],
            
            os.system(wgsmetrics_cmd)
        else:
            hsmetrics_cmd = ("%s -jar %s CalculateHsMetrics BI=%s TI=%s" 
                         "METRIC_ACCUMULATION_LEVEL=ALL_READS I=%s O=%s"
                         "MQ=20 Q=10"%(self.args['java_location'],
                                       self.args['picard_location'],
                                       self.args['bait_file'],
                                       self.args['target_file'],
                                       self.args['bam_file'],
                                       self.output_file))
            ## Run the command
            os.system(hsmetrics_cmd)

    
class RootTask(luigi.WrapperTask):
    """
    Root Task for the pipeline. Entry point for
    dependency graph 
    """
    
    args = luigi.DictParameter() ## User defined commandline args
    
    
    def requires(self):
        return [RunMetrics()]
