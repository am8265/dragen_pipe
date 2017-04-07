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

class WrapperEthnicity(luigi.WrapperTask):
    """ Update the database with ethnicity probabilities
    """

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    workdir = luigi.Parameter()

    def __init__(self,*args,**kwargs):
        super(WrapperEthnicity,self).__init__(*args,**kwargs)
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.log_file = os.path.join(self.workdir,"update_db_ethnicity.log")
        logging.basicConfig(filename=self.log_file,level=logging.DEBUG)
        
    def requires(self):
        for chgvid,pseudo_prepid,seqtype,alignloc in get_samples_to_predict():
            yield PredictAndUpdate(sample_name=chgvid,pseudo_prepid=pseudo_prepid,
                                   sample_type=seqtype,alignseqfileloc=alignloc,
                                   workdir = self.workdir,run_locally=True)

class PredictAndUpdate(SGEJobTask):
    """ Run the classifier and get the probablities
    """
    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"
    
    sample_name = luigi.Parameter()
    pseudo_prepid = luigi.Parameter()
    sample_type = luigi.Parameter()
    alignseqfileloc = luigi.Parameter()
    workdir = luigi.Parameter()
    training_model = luigi.Parameter(default="{0}/data/37MB_markers_model.obj".format(os.getcwd()))
    training_snps = luigi.Parameter(default="{0}/data/filtered.37MB.master.training.map".format(os.getcwd()))
    ethnicity_wrapper = luigi.Parameter(default="{0}/model_ancestry.py".format(os.getcwd()))
    
    def __init__(self,*args,**kwargs):
        """
        """
        super(PredictAndUpdate,self).__init__(*args,**kwargs)
        self.outdir = os.path.join(self.alignseqfileloc,'{0}.{1}'.format(self.sample_name,self.pseudo_prepid))
        self.output_probs = os.path.join(self.outdir,"ethnicity.probs.txt")
        self.input_ped = os.path.join(self.outdir,self.sample_name+'.'+str(self.pseudo_prepid)+'.ped')        
        
    def work(self):
        """
        Run the python script to predict ethnicities
        """
        cmd = ( """ python {0} --output-prob-file {1}"""
                """ --testing-ped {2} --mapfile {3}"""
                """ --input-model {4} """.format(
                    self.ethnicity_wrapper, self.output_probs,self.input_ped,
                    self.training_snps, self.training_model)
                )        
        run_shellcmd(cmd)
        update_probs(self.output_probs,self.pseudo_prepid)
        update_ped_status(self.pseudo_prepid,"predict")
        
    def requires(self):
        """
        """        
        return CreatePed(pseudo_prepid = self.pseudo_prepid,sample_name=self.sample_name,
                         alignseqfile=self.alignseqfileloc,
                         sample_type=self.sample_type,markers=self.training_snps,
                         run_locally=True)
    
    def output(self):
        """
        """        
        return SQLTarget(self.pseudo_prepid,"predict")       

