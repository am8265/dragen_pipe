#!/usr/bin/env Python
import sys
import subprocess
import os
import utils 
import luigi


class CreateTargetFile(luigi.Task):
    """ 
    A Luigi Task
    """

    ## Parameters for this task 
    bed_file = luigi.Parameter()
    java_location = luigi.Parameter()
    picard_location = luigi.Parameter()
    seqdict_file = luigi.Parameter()
    output_dir = luigi.Parameter()
    
    ## Prepare the output name for the TargetFile
    output_name = utils.get_filename_from_fullpath(self.bed_file)
    output_file = self.output_dir + self.output_name + ".list"
    
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

        bedtointerval_cmd = ("%s -jar %s BedToIntervalList I=%s SD=%s 
                             OUTPUT=%s"%(self.java_location,
                                         self.picard_location,
                                         self.bed_file,
                                         self.seqdict_file,
                                         self.output_file))
        os.system(bedtointerval_cmd)

class RunHsMetrics(luigi.Task):
    """ 
    A luigi task
    """

    ## Parameters for this task
    java_location = luigi.Parameter()
    picard_location = luigi.Parameter()
    bait_file = luigi.Parameter()
    target_file = luigi.Parameter()
    bam_file = luigi.Parameter()
    sample = luigi.Parameter()
    output_dir = luigi.Parameter()
    create_targetfile = luigi.BoolParameter()
    
    ## Prepare the outputfile name
    output_file = self.output_dir + self.sample + ".hsmetrics.txt"

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
        if self.create_targetfile == True:
            return CreateTargetFile()
    
    def run(self):
        """
        Run Picard CalculateHsMetrics
        """
        
        hsmetrics_cmd = ("%s -jar %s CalculateHsMetrics BI=%s TI=%s 
                         METRIC_ACCUMULATION_LEVEL=ALL_READS I=%s O=%s
                         MQ=20 Q=10"%(java_location,picard_location,
                                      bait_file,target_file,bam_file,
                                      output_file))
        ## Run the command
        os.system(hsmetrics_cmd)

        
