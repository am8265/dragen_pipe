### Python specific modules

import MySQLdb
from automated_relatedness import get_connection
import os
import subprocess
import luigi
from luigi.contrib.sge import SGEJobTask
from datetime import datetime
import logging
from luigi.task import Task
### Modules from other scripts in this directory
# from automated_ethnicity_dragen import get_samples_to_predict, update_probs
from automated_relatedness import run_shellcmd, update_ped_status
from create_ped_luigi_dragen import CreatePed, SQLTarget

base_directory = os.path.dirname(os.path.abspath(__file__))
###################### start of automated_ethnicity_dragen.py

def update_probs(output_prob_file,pseudo_prepid):
    db = get_connection(db="seqdb")
    cmd = ( """ UPDATE dragen_qc_metrics SET Caucasian_prob = {0},MiddleEastern_prob = {1},"""
      """ Hispanic_prob = {2}, EastAsian_prob = {3} , SouthAsian_prob = {4},"""
      """ African_prob = {5}, genotyping_rate = {6} WHERE pseudo_prepid = {7} """ )
    
    with open(output_prob_file,'r') as IN:
        for line in IN:
            line = line.strip('\n')
            caucasian,middleeast,hispanic,eastasian,southasian,african,genorate = line.split('\t')
            cur = db.cursor()
            db_cmd = cmd.format(caucasian,middleeast,hispanic,eastasian,southasian,african,genorate,pseudo_prepid)
            print db_cmd
            cur.execute(db_cmd)
            db.commit()
        db.close()
        
def get_samples_to_predict():
    ######## CLEARLY, do NOT use the silly ped table!?!
    ######## CLEARLY, do NOT use the silly ped table!?!
    ######## CLEARLY, do NOT use the silly ped table!?!
    max_samples = 1
    query = "select sample_name, d.pseudo_prepid, sample_type, alignseqfileloc from dragen_sample_metadata d join dragen_qc_metrics q on d.pseudo_prepid=q.pseudo_prepid \
      where genotyping_rate is null and is_merged in (40,100) order by pseudo_prepid desc limit {}".format(max_samples);
    # print("\n\n{}\n\n".format(query))
    # query = "select sample_name, d.pseudo_prepid, sample_type, is_merged, alignseqfileloc, genotyping_rate from dragen_sample_metadata d join dragen_qc_metrics q on d.pseudo_prepid=q.pseudo_prepid where genotyping_rate is null and is_merged in (40,100) order by pseudo_prepid desc limit {}".format(max_samples);
    db = get_connection(db="seqdb")
    cur = db.cursor()
    cur.execute(query)
    result = cur.fetchall()
    ##### dragen_ped_status is JUST a pipeline table for each of the stages and separate cols!?!?!?!?!?!?
    ##### presumably populates in one go and then tries to run them in one batch and if anything happens it all breaks.
    ##### may need to wipe --base-output-directory and truncate table!?!?
    # max_samples = 8000
    i = 0 
    for res in result[0:max_samples]:
        i+=1
        yield res
  

###################### end of automated_ethnicity_dragen.py

# class PredictAndUpdate(Task):
class PredictAndUpdate(SGEJobTask):

    # Elements outside the __init__ method are static elements
    base_output_directory = luigi.OutputDirectoryParameter()
    # output_directory = luigi.OutputDirectoryParameter()
    pseudo_prepid = luigi.IntParameter()


    training_model      = luigi.InputFileParameter( default=os.path.join( base_directory, "data", "37MB_markers_model.obj") )
    markers             = luigi.InputFileParameter( default=os.path.join( base_directory, "data", "filtered.37MB.master.training.map") )
    ethnicity_wrapper   = luigi.InputFileParameter( default=os.path.join( base_directory, "model_ancestry.py") )
    
    # Elements inside the __init__ method are elements of the object (self)
    def __init__(self,*args,**kwargs):

        super(PredictAndUpdate,self).__init__(*args,**kwargs)

        query="select sample_name, d.pseudo_prepid, sample_type, alignseqfileloc from dragen_sample_metadata d join dragen_qc_metrics q on d.pseudo_prepid=q.pseudo_prepid \
          where genotyping_rate is null and is_merged in (40,100) "
        if self.pseudo_prepid!=0 :
            query += "and d.pseudo_prepid = {}".format(self.pseudo_prepid)
        else:
            query += "order by pseudo_prepid desc limit 1"

        print("using {}".format(query))

        db = get_connection(db="seqdb")
        cur = db.cursor()
        cur.execute(query)
        if cur.rowcount==0:
            print("there are no matching samples")
            exit(1)

        # self.run_locally=True
        # from pprint import pprint as pp
        # pp(vars(self))

        self.sample_name, self.pseudo_prepid, self.sample_type, align_loc = cur.fetchone()

        self.align_loc  =   os.path.join( align_loc, "{sample_name}.{pseudo_prepid}".format( sample_name=self.sample_name, pseudo_prepid=self.pseudo_prepid ) )

        self.output_directory=os.path.join( self.base_output_directory, "{sample_name}.{pseudo_prepid}.{sample_type}".format( 
          sample_name=self.sample_name, pseudo_prepid=self.pseudo_prepid, sample_type=self.sample_type ) ) 
        # sample_name = luigi.Parameter()
        # pseudo_prepid = luigi.IntParameter()
        # sample_type = luigi.Parameter()
        # align_loc = luigi.InputDirectoryParameter()

        self.run_locally=True

        if not os.path.isdir(self.output_directory):
            os.mkdir(self.output_directory)
        else:
            print("warning: dir '{}' already exists".format(self.output_directory))

        cur.execute("update dragen_qc_metrics set genotyping_rate = -1.0 where genotyping_rate is null and pseudo_prepid = {}".format(self.pseudo_prepid))
        print("we modified {} rows".format(cur.rowcount))
        if cur.rowcount != 1:
            print("unable to get lock on pseudo_prepid = {}".format(self.pseudo_prepid))
        db.commit()
        self.output_probs = os.path.join( self.output_directory ,"ethnicity.probs.txt")
        self.input_ped = os.path.join( self.output_directory, "{sample_name}.{pseudo_prepid}.ped".format(
          sample_name=self.sample_name, pseudo_prepid=self.pseudo_prepid ) )
 
    def work(self):
        """Run the python script to predict ethnicities"""
        cmd = ( " /nfs/goldstein/software/python2.7.7/bin/python {0} --output-prob-file {1} --testing-ped {2} --mapfile {3} --input-model {4} """.format(
          self.ethnicity_wrapper, self.output_probs,self.input_ped, self.markers, self.training_model ) )        
        print("ethnicity_check_dragen.py '{}'".format(cmd))
        run_shellcmd(cmd)

        update_probs(self.output_probs, self.pseudo_prepid)

        update_ped_status(self.pseudo_prepid, "predict")
        
    def requires(self):
        # return self.clone(CreatePed)
        return CreatePed(sample_name=self.sample_name, pseudo_prepid=self.pseudo_prepid, sample_type=self.sample_type, 
          align_loc=self.align_loc, output_directory=self.output_directory, markers=self.markers)
    
    def output(self):
        return SQLTarget(self.pseudo_prepid,"predict")       
