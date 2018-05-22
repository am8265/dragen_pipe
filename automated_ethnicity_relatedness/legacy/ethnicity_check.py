import MySQLdb
import os
import subprocess
import luigi
from luigi.contrib.sge import SGEJobTask
from datetime import datetime
from automated_ethnicity import update_ped_status,get_samples_to_predict, update_probs
from automated_relatedness import MyExtTask,get_connection,get_familyid,run_shellcmd
import logging

base_directory = os.path.dirname(os.path.abspath(__file__))

class WrapperEthnicity(luigi.WrapperTask):
    """ Update the database with ethnicity probabilities
    """

    base_output_directory = luigi.OutputDirectoryParameter()

    def __init__(self,*args,**kwargs):
        super(WrapperEthnicity,self).__init__(*args,**kwargs)
        
    def requires(self):
        # for f' sake...?!?
        for sample_name, prepid, sample_type, align_loc, priority in get_samples_to_predict():
            yield PredictAndUpdate(
                sample_name=sample_name, prepid=prepid, sample_type=sample_type,
                align_loc=align_loc, output_directory=os.path.join(
                    self.base_output_directory, "{sample_name}.{prepid}.{sample_type}".
                    format(sample_name=sample_name, prepid=prepid, sample_type=sample_type)),
                bioinfo_priority=int(priority))

class PredictAndUpdate(SGEJobTask):
    """
    Run the classifier and get the probablities
    """

    sample_name = luigi.Parameter()
    prepid = luigi.IntParameter()
    sample_type = luigi.Parameter()
    align_loc = luigi.InputDirectoryParameter()
    output_directory = luigi.OutputDirectoryParameter()
    bioinfo_priority = luigi.IntParameter()
    training_model = luigi.InputFileParameter(
        default="/nfs/seqscratch_ssd/bc2675/annodb_ethnicities/"
        "37MB_markers_model.obj")
    markers = luigi.InputFileParameter(
        default="/nfs/seqscratch_ssd/bc2675/annodb_ethnicities/"
        "filtered.37MB.master.training.map")
    
    def __init__(self,*args,**kwargs):
        super(PredictAndUpdate,self).__init__(*args,**kwargs)
        self.output_probs = os.path.join(
            self.output_directory,"ethnicity.probs.txt")
        self.input_ped = os.path.join(
            self.output_directory, "{sample_name}.{prepid}.ped".format(
                sample_name=self.sample_name, prepid=self.prepid))

    @property
    def priority(self):
        return 10 - self.bioinfo_priority
        
    def work(self):
        """ Run the python script to predict ethnicities
        """
        cmd = ( """ python {0} --output-prob-file {1}"""
                """ --testing-ped {2} --mapfile {3}"""
                """ --input-model {4} """.format(
                    self.ethnicity_wrapper, self.output_probs,self.input_ped,
                    self.markers, self.training_model)
                )        
        run_shellcmd(cmd)
        update_probs(self.output_probs,self.prepid,
                     self.sample_name,self.sample_type)
        update_ped_status(self.sample_name,self.prepid,
                          self.sample_type,"predict")
        
    def requires(self):
        """
        """        
        return self.clone(CreatePed)
    
    def output(self):
        """
        """        
        return SQLTarget(self.prepid,self.sample_name,
                         self.sample_type,"predict")
        

class CreatePed(SGEJobTask):
    """ Create Ped 
    """

    output_directory = luigi.OutputDirectoryParameter()
    prepid = luigi.IntParameter()
    sample_name = luigi.Parameter()
    align_loc = luigi.InputDirectoryParameter()
    sample_type = luigi.Parameter()
    markers = luigi.InputFileParameter()
    bioinfo_priority = luigi.IntParameter()
    
    def __init__(self,*args,**kwargs):
        """ Initialize the class
        """
        super(CreatePed,self).__init__(*args,**kwargs)
        archive_loc = os.path.join(self.align_loc,'combined')
        self.bam = os.path.join(
            archive_loc, "combined_rmdup_realn_recal.bam")
        vcf = os.path.join(
            archive_loc, self.sample_name + ".analysisReady.vcf")
        if not os.path.isfile(vcf):
            vcf = vcf + ".gz"
        if not os.path.isfile(vcf):
            vcf = os.path.join(
                archive_loc, self.sample_name + ".analysisReady.annotated.vcf")
        if not os.path.isfile(vcf):
            vcf = vcf + ".gz"
        if not os.path.isfile(vcf):
            raise OSError("Couldn't find VCF at {loc}".format(loc=archive_loc))
        self.vcf = vcf
        self.log_file = os.path.join(self.output_directory,"{0}.{1}.createped.log".format(self.sample_name,self.prepid))
        self.famid = get_familyid(self.sample_name)
        self.pedmap_script = "/nfs/seqscratch_ssd/bc2675/annodb_ethnicities/create_ped_map.py"
        self.cmd="/nfs/goldstein/software/python2.7.7/bin/python {0} --vcf {1} --sample_id {2} --markers {3} --seqtype {4} --pseudo_prepid {5} --bam {6} --stem {2}.{5} -out {7} --logfile {8} --family_id {9}".format(
          self.pedmap_script,self.vcf,self.sample_name,self.markers,self.sample_type,self.prepid,self.bam,self.output_directory,self.log_file,self.famid
        )
        
    @property
    def priority(self):
        return 10 - self.bioinfo_priority

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
        return [MyExtTask(self.bam),MyExtTask(self.vcf)]

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
