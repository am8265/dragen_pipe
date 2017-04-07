import os
import luigi
from luigi.contrib.sge import SGEJobTask
from automated_relatedness import MyExtTask, get_connection, update_ped_status, run_shellcmd, get_familyid

class CreatePed(SGEJobTask):
    """ Create Ped 
    """

    n_cpu = 1
    parallel_env =  "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    pseudo_prepid = luigi.Parameter()
    sample_name = luigi.Parameter()
    alignseqfile = luigi.Parameter()
    sample_type = luigi.Parameter()
    markers = luigi.Parameter()
    
    def __init__(self,*args,**kwargs):
        """ Initialize the task class
        """
        super(CreatePed,self).__init__(*args,**kwargs)
        self.outdir = os.path.join(self.alignseqfile,'{0}.{1}'.format(self.sample_name,self.pseudo_prepid))
        self.bam = os.path.join(self.outdir,'{0}.{1}.realn.recal.bam'.format(self.sample_name,self.pseudo_prepid))
        self.vcf_gz = os.path.join(self.outdir,'{0}.{1}.analysisReady.annotated.vcf.gz'.format(self.sample_name,self.pseudo_prepid))
        self.log_file = os.path.join(self.outdir,'logs/{0}.{1}.createped.log'.format(self.sample_name,self.pseudo_prepid))
        self.famid = get_familyid(self.sample_name)
        self.pedmap_script = os.path.join(os.getcwd(),"create_ped_map.py")
        self.cmd="/nfs/goldstein/software/python2.7.7/bin/python {0} --vcf {1} --sample_id {2} --markers {3} --seqtype {4} --pseudo_prepid {5} --bam {6} --stem {2}.{5} -out {7} --logfile {8} --family_id {9}".format(self.pedmap_script,self.vcf_gz,self.sample_name,self.markers,self.sample_type,self.pseudo_prepid,self.bam,self.outdir,self.log_file,self.famid)
        
    def work(self):
        """ The command(s) to run 
        """        
        run_shellcmd(self.cmd)
        ## Update the database with correct status
        update_ped_status(self.pseudo_prepid,"create_ped")
        
    def output(self):
        """ The output from this task 
        """        
        return SQLTarget(self.pseudo_prepid,"create_ped")
        
    def requires(self):
        """ THe dependencies 
        """        
        return [MyExtTask(self.bam),MyExtTask(self.vcf_gz)]

class SQLTarget(luigi.Target):
    """ A luigi target class describing verification of the entries in the database
    """    
    def __init__(self, pseudo_prepid, field):
        self.pseudo_prepid = pseudo_prepid
        self.field = field
        
    def exists(self):
        db = get_connection(db="seqdb")
        try:
            cmd = (""" SELECT {0} FROM dragen_ped_status """
                   """ WHERE pseudo_prepid = {1}""".format(self.field,self.pseudo_prepid)
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
                cmd = (""" INSERT INTO dragen_ped_status"""
                       """ (pseudo_prepid)"""
                       """ VALUES"""
                       """ ({0});""".format(self.pseudo_prepid)
                       )
                cur = db.cursor()
                cur.execute(cmd)
                db.commit()
                db.close()
                return False 
        finally:
            if db.open:
                db.close()
