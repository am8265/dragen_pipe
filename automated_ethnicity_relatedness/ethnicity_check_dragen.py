### Python specific modules
import MySQLdb
import os
import subprocess
import luigi
from luigi.contrib.sge import SGEJobTask
from datetime import datetime
import logging
### Modules from other scripts in this directory
from automated_ethnicity_dragen import get_samples_to_predict, update_probs
from automated_relatedness import MyExtTask,get_connection,get_familyid,run_shellcmd, update_ped_status
from create_ped_luigi_dragen import CreatePed, SQLTarget

base_directory = os.path.dirname(os.path.abspath(__file__))

class WrapperEthnicity(luigi.WrapperTask):
    """ Update the database with ethnicity probabilities
    """

    base_output_directory = luigi.OutputDirectoryParameter()

    def __init__(self,*args,**kwargs):
        super(WrapperEthnicity,self).__init__(*args,**kwargs)
        
    def requires(self):
        # for f' sake?!?
        for sample_name, pseudo_prepid, sample_type, align_loc in get_samples_to_predict():
            yield self.clone(PredictAndUpdate, sample_name=sample_name,
                             pseudo_prepid=pseudo_prepid,
                             sample_type=sample_type,
                             align_loc=os.path.join(align_loc, "{sample_name}.{pseudo_prepid}".format(
                                 sample_name=sample_name, pseudo_prepid=pseudo_prepid)),
                             output_directory=os.path.join(
                self.base_output_directory, "{sample_name}.{pseudo_prepid}.{sample_type}".format(
                    sample_name=sample_name, pseudo_prepid=pseudo_prepid,
                    sample_type=sample_type)))

class PredictAndUpdate(SGEJobTask):
    """ Run the classifier and get the probablities
    """
    
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.IntParameter()
    sample_type = luigi.Parameter()
    align_loc = luigi.InputDirectoryParameter()
    output_directory = luigi.OutputDirectoryParameter()
    training_model = luigi.InputFileParameter(default=os.path.join(
        base_directory, "data", "37MB_markers_model.obj"))
    markers = luigi.InputFileParameter(default=os.path.join(
        base_directory, "data", "filtered.37MB.master.training.map"))
    ethnicity_wrapper = luigi.InputFileParameter(default=os.path.join(
        base_directory, "model_ancestry.py"))
    
    def __init__(self,*args,**kwargs):
        super(PredictAndUpdate,self).__init__(*args,**kwargs)
        self.run_locally=True
        os.mkdir(self.output_directory)
        self.output_probs = os.path.join(
            self.output_directory ,"ethnicity.probs.txt")
        self.input_ped = os.path.join(
            self.output_directory, "{sample_name}.{pseudo_prepid}.ped".format(
                sample_name=self.sample_name, pseudo_prepid=self.pseudo_prepid))
        
    def work(self):
        """
        Run the python script to predict ethnicities
        """
        cmd = ( """ /nfs/goldstein/software/python2.7.7/bin/python {0} --output-prob-file {1}"""
                """ --testing-ped {2} --mapfile {3}"""
                """ --input-model {4} """.format(
                    self.ethnicity_wrapper, self.output_probs,self.input_ped,
                    self.markers, self.training_model)
                )        
        run_shellcmd(cmd)
        update_probs(self.output_probs, self.pseudo_prepid)
        update_ped_status(self.pseudo_prepid, "predict")
        
    def requires(self):
        """
        """        
        return self.clone(CreatePed)
    
    def output(self):
        """
        """        
        return SQLTarget(self.pseudo_prepid,"predict")       

