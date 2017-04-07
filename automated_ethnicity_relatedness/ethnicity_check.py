import MySQLdb
import os
import subprocess
import luigi
from luigi.contrib.sge import SGEJobTask
from datetime import datetime
from automated_ethnicity import update_ped_status,get_samples_to_predict, update_probs
from automated_relatedness import MyExtTask,get_connection,get_familyid,run_shellcmd
import logging


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
        for chgvid,prepid,seqtype,alignloc in get_samples_to_predict():
            yield PredictAndUpdate(sample_name=chgvid,prepid=prepid,
                                   sample_type=seqtype,alignseqfileloc=alignloc,
                                   workdir = self.workdir,run_locally=True)

class PredictAndUpdate(SGEJobTask):
    """
    Run the classifier and get the probablities
    """

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    sample_name = luigi.Parameter()
    prepid = luigi.Parameter()
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
        self.outdir = '/nfs/goldsteindata/ethnicity_predictions/{0}.{1}.{2}'.format(
            self.sample_name,self.prepid,self.sample_type)
        self.output_probs = os.path.join(self.outdir,"ethnicity.probs.txt")
        self.input_ped = os.path.join(self.outdir,self.sample_name+'.'+str(self.prepid)+'.ped')
        self.testing_geno_file = os.path.join(self.workdir,self.sample_name+'.'+str(self.prepid)+'.'+self.sample_type+'.geno')
        
    def work(self):
        """ Run the python script to predict ethnicities
        """
        cmd = ( """ python {0} --output-prob-file {1}"""
                """ --testing-ped {2} --mapfile {3}"""
                """ --input-model {4} """.format(
                    self.ethnicity_wrapper, self.output_probs,self.input_ped,
                    self.training_snps, self.training_model)
                )        
        run_shellcmd(cmd)
        update_probs(self.output_probs,self.prepid,
                     self.sample_name,self.sample_type)
        update_ped_status(self.sample_name,self.prepid,
                          self.sample_type,"predict")
        
    def requires(self):
        """
        """        
        return CreatePed(prepid = self.prepid,sample_name=self.sample_name,
                         alignseqfile=self.alignseqfileloc,
                         sample_type=self.sample_type,
                         markers=self.training_snps,run_locally=True)
    
    def output(self):
        """
        """        
        return SQLTarget(self.prepid,self.sample_name,
                         self.sample_type,"predict")
        

class CreatePed(SGEJobTask):
    """ Create Ped 
    """

    n_cpu = 1
    parallel_env =  "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    prepid = luigi.Parameter()
    sample_name = luigi.Parameter()
    alignseqfile = luigi.Parameter()
    sample_type = luigi.Parameter()
    markers = luigi.Parameter()
    
    def __init__(self,*args,**kwargs):
        """ Initialize the class
        """
        super(CreatePed,self).__init__(*args,**kwargs)
        self.outdir = '/nfs/goldsteindata/ethnicity_predictions/{0}.{1}.{2}'.format(self.sample_name,self.prepid,self.sample_type)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        self.archive_loc = os.path.join(self.alignseqfile,'combined')
        self.bam = os.path.join(self.archive_loc,self.sample_name+'_final.bam')
        self.vcf_gz = os.path.join(self.archive_loc,self.sample_name+'.analysisReady.annotated.vcf.gz')
        self.log_file = os.path.join(self.outdir,"{0}.{1}.createped.log".format(self.sample_name,self.prepid))
        self.famid = get_familyid(self.sample_name)
        self.pedmap_script = os.path.join(os.getcwd(),"create_ped_map.py")
        self.cmd="/nfs/goldstein/software/python2.7.7/bin/python {0} --vcf {1} --sample_id {2} --markers {3} --seqtype {4} --pseudo_prepid {5} --bam {6} --stem {2}.{5} -out {7} --logfile {8} --family_id {9}".format(self.pedmap_script,self.vcf_gz,self.sample_name,self.markers,self.sample_type,self.prepid,self.bam,self.outdir,self.log_file,self.famid)
        
    def work(self):
        """
        """        
        run_shellcmd(self.cmd)
        ## Update the database with correct status
        update_ped_status(self.sample_name,self.prepid,self.sample_type,"create_ped")
        
    def output(self):
        """
        """        
        return SQLTarget(self.prepid,self.sample_name,self.sample_type,"create_ped")
        
    def requires(self):
        """
        """        
        return [MyExtTask(self.bam),MyExtTask(self.vcf_gz)]

class SQLTarget(luigi.Target):
    """ A luigi target class describing verification of the entries in the database
    """    
    def __init__(self, prepid,chgvid,seqtype,field):
        self.prepid = prepid
        self.chgvid = chgvid
        self.seqtype = seqtype
        self.field = field
        
    def exists(self):
        db = get_connection(db="seqdb")
        try:
            cmd = (""" SELECT {0} FROM ethnicity_status """
                   """ WHERE prepid = {1} AND CHGVID = '{2}' """
                   """ AND SeqType = '{3}' """.format(self.field,self.prepid,
                                                    self.chgvid,self.seqtype)
                   )
            cur = db.cursor()
            cur.execute(cmd)
            row = cur.fetchone()
            if row:
                if row[0] == 1:
                    return True
                else:
                    return False
            else: ## Not yet initialized in the status table, so initialize it here
                cmd = (""" INSERT INTO ethnicity_status"""
                       """ (prepid,CHGVID,SeqType)"""
                       """ VALUES"""
                       """ ({0},'{1}','{2}');""".format(self.prepid,self.chgvid,
                                                    self.seqtype))
                cur = db.cursor()
                cur.execute(cmd)
                db.commit()
                db.close()
                return False 
        finally:
            if db.open:
                db.close()
