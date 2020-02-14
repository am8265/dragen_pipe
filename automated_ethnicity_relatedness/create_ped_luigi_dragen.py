import os
import luigi
from luigi.contrib.sge import SGEJobTask
from automated_relatedness import MyExtTask, get_connection, update_ped_status, run_shellcmd, get_familyid

base_directory = os.path.dirname(os.path.abspath(__file__))

class CreatePed(SGEJobTask):

    pseudo_prepid = luigi.IntParameter()
    sample_name = luigi.Parameter()
    align_loc = luigi.InputDirectoryParameter()
    output_directory = luigi.OutputDirectoryParameter()
    sample_type = luigi.Parameter()
    markers = luigi.InputFileParameter()
    
    def __init__(self,*args,**kwargs):
        super(CreatePed,self).__init__(*args,**kwargs)
        self.run_locally=True
        self.bam = os.path.join(self.align_loc,'{0}.{1}.realn.recal.bam'.format(self.sample_name,self.pseudo_prepid))
        self.vcf_gz = os.path.join(self.align_loc,'{0}.{1}.analysisReady.annotated.vcf.gz'.format(self.sample_name,self.pseudo_prepid))
        self.log_file = os.path.join(self.output_directory,'{0}.{1}.createped.log'.format(self.sample_name,self.pseudo_prepid))
        # self.log_file = os.path.join(self.align_loc,'logs/{0}.{1}.createped.log'.format(self.sample_name,self.pseudo_prepid))
        self.famid = get_familyid(self.sample_name)
        self.pedmap_script = os.path.join(base_directory, "create_ped_map.py")
        self.cmd=(  "/nfs/goldstein/software/python2.7.7/bin/python {0} --vcf {1} --sample_id {2} --markers {3} --seqtype {4} --pseudo_prepid {5} --bam "
                    "{6} --stem {2}.{5} -out {7} --logfile {8} --family_id {9}".format( self.pedmap_script, self.vcf_gz, self.sample_name, self.markers,
                       self.sample_type, self.pseudo_prepid, self.bam, self.output_directory, self.log_file, self.famid ) )
        self.final_ped = '{0}/{1}.{2}.ped'.format(self.output_directory,self.sample_name,self.pseudo_prepid,self.sample_type)
        print("\n\nget cmd from init\n{}\noutput\n{}\n\n".format(self.cmd,self.final_ped))
        db = get_connection(db="seqdb")
        cmd = ("SELECT create_ped FROM dragen_ped_status WHERE pseudo_prepid = {}".format(self.pseudo_prepid) )
        cur = db.cursor()
        cur.execute(cmd)
        row = cur.fetchone()
        if row == None:
            print("we need to insert it for now");
            cmd = (" INSERT INTO dragen_ped_status (pseudo_prepid) VALUES ({0});".format(self.pseudo_prepid) )
            cur = db.cursor()
            print("running '{}'".format(cmd))
            cur.execute(cmd)
            db.commit()
            db.close()
            if os.path.exists(self.final_ped):
                print("there's already a ped file {}".format(self.final_ped))
                update_ped_status(self.pseudo_prepid,"create_ped")
            else:
                print("there's no ped file {}".format(self.final_ped))
        else:
            print("already inserted")
        # SQLTarget(self.pseudo_prepid,"create_ped").exists()
        
    def work(self):
        db = get_connection(db="seqdb")
        try:
            cur = db.cursor()
            cmd = ("select * from dragen_ped_status where pseudo_prepid = {}".format(self.pseudo_prepid) )
            print("locking")
            cur.execute(cmd)
            from pprint import pprint as pp
            pp(cur.fetchall())
            cmd = ("update dragen_ped_status set create_ped = -1 WHERE create_ped = 0 and pseudo_prepid = {}".format(self.pseudo_prepid) )
            print("locking")
            cur.execute(cmd)
            print("got '{}'".format(cur.rowcount))
            if cur.rowcount != 1 :
                print("unable to lock")
                exit(1);
            db.commit()
            db.close()
        except:
            print("unable to lock")
            exit(1);
        print("running '{}'".format(cmd))
        print("\n\ncreate_ped_luigi_dragen.py '{}'\n\n".format(self.cmd))
        run_shellcmd(self.cmd)
        ## Update the database with correct status
        update_ped_status(self.pseudo_prepid,"create_ped")
        
    def output(self):
        # return luigi.LocalTarget(self.final_ped)
        return [luigi.LocalTarget(self.final_ped),SQLTarget(self.pseudo_prepid,"create_ped")]
        # return SQLTarget(self.pseudo_prepid,"create_ped")
        
    def requires(self):
        return [MyExtTask(self.bam),MyExtTask(self.vcf_gz)]

class SQLTarget(luigi.Target):
    def __init__(self, pseudo_prepid, field):
        print(">>> silly sqltarget init")
        self.pseudo_prepid = pseudo_prepid
        self.field = field
        
    def exists(self):
        print("  >>> silly sqltarget exists")
        db = get_connection(db="seqdb")
        try:
            cmd = ("SELECT {0} FROM dragen_ped_status WHERE pseudo_prepid = {1}".format(self.field,self.pseudo_prepid) )
            cur = db.cursor()
            # print("running '{}'".format(cmd))
            cur.execute(cmd)
            row = cur.fetchone()
            from pprint import pprint as pp
            # pp(row)
            if row:
                if row[0] == 1:
                    print("     >>> sqltarget already done\n")
                    return True
                else:
                    return False
            else: ## Not yet initialized in the status table, so initialize it here
                cmd = (" INSERT INTO dragen_ped_status (pseudo_prepid) VALUES ({0});".format(self.pseudo_prepid) )
                cur = db.cursor()
                # print("running '{}'".format(cmd))
                cur.execute(cmd)
                db.commit()
                db.close()
                return False 
        finally:
            if db.open:
                db.close()
