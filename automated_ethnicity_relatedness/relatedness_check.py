### Python specific modules
import MySQLdb
import os
from ConfigParser import ConfigParser
import subprocess
import sys
import luigi
from luigi.contrib.sge import SGEJobTask
from datetime import datetime
### Modules from other scripts in this directory
from automated_relatedness import get_connection, get_samples_append, update_ped_status, run_shellcmd
from create_ped_luigi_dragen import CreatePed, SQLTarget

###### UpdateDBRelatedness & RunKing use silly brute-force (not accomodating sample list...) touch files while AppendMasterPed uses the append_ped field
###### question is how do we go about running this stuff - all of it, by batch, all newer samples?!?
###### perhaps just load the .ped data?!?
######
###### cp /nfs/goldstein/software/dragen_pipe/master/automated_ethnicity_relatedness/data/filtered.37MB.master.training.map /nfs/fastq_temp/dsth/relatedness/combined.map
###### :> /nfs/fastq_temp/dsth/relatedness/combined.ped 
###### rm /nfs/fastq_temp/dsth/relatedness/*.success
###### mysql_seqdb_q "update dragen_ped_status set append_ped = 0"
######
###### need to re-set create_ped where appropriate!?!?
######
###### 'should' run on ALL non-custom_capture samples?!?
######

class UpdateDBRelatedness(SGEJobTask):

    # n_cpu           = 1
    # parallel_env    = "threaded"
    # shared_tmp_dir  = "/nfs/seqscratch09/tmp/luigi_test"

    batch_size = luigi.IntParameter()
    ######## the f'ing map?!
    masterped       = luigi.Parameter(default="/nfs/fastq_temp/dsth/relatedness/combined.ped")
    output_directory= luigi.Parameter()
    markers         = luigi.Parameter(default="{0}/filtered.37MB.master.training.map".format(os.getcwd()))
    genotyping_dir  = luigi.Parameter()

    def __init__(self,*args,**kwargs):

        super(UpdateDBRelatedness,self).__init__(*args,**kwargs)

        self.log_file           = os.path.join( self.output_directory,"update_db_relatedness.log" )
        self.temp_bed           = os.path.join( self.output_directory,"temp_masterped" )
        self.output_prefix      = os.path.join( self.output_directory,"relatedness" )
        self.king_btw           = self.output_prefix+".kin0"
        self.king_wtn           = self.output_prefix+".kin"
        self.verification_file  = os.path.join( self.output_directory,"update_db.success" )
        self.btw_fam_threshold  = 0.1
        self.run_locally=True
        
    def work(self):
        self.update_between_family()
        self.update_within_family()
        os.system("touch {0}".format(self.verification_file))
        
    def requires(self):
        return self.clone(RunKing)
        # return RunKing( masterped=self.masterped,output_directory=self.output_directory, markers=self.markers, run_locally=True, genotyping_dir=self.genotyping_dir )
        
    def output(self):
        #return MyExtTask(self.verification_file)
        return luigi.LocalTarget(self.verification_file)

    def update_between_family(self):
        """ Check and update relatedness table
        between families/unrealated individuals
        """

        statement = ("""INSERT INTO between_fam_relatedness """
                     """(sample_name1, pseudo_prepid1, FamilyID1 """
                     """,sample_name2, pseudo_prepid2, FamilyID2, coef) """
                     """ VALUES ('{0}',{1},'{2}','{3}',{4},'{5}',{6}) ON DUPLICATE KEY """
                     """ UPDATE coef = {6} """)
        
        if self.parse_kin_btw(self.king_btw): ## The file was not empty , an empty file means no relatedness was inferred in
            ## individuals between families 
            for sample_name,pseudo,family_name,related_sample_name,related_pseudo,related_family_name,coef in self.parse_kin_btw(self.king_btw):
                self.update_db(statement.format(sample_name,pseudo,family_name,related_sample_name,
                                           related_pseudo,related_family_name,coef))


    def check_relation(self,relation,coef):
        """
        Implments relation to proband guidlines from igm bioinformatics wiki
        relation : str ; The relation to proband
        coef : float ; the kinship coefficient 
        """

        if relation in ["Great-great grandparent","Great-great grandchild",
                        "Second cousin","Third cousin","Fourth cousin",
                        "First cousin once removed","First cousin twice removed",
                        "Second cousin once removed","Second cousin twice removed",
                        "Thrid cousin once removed","Half first cousin",
                        "Spouse","Step sibling","Step parent","Step child",
                        "Adoptive family member"]:
            return coef <= 0.1
        
        elif relation in ["Great grandparent","Great grandchild","Great aunt-uncle",
                          "Great niece-nephew","Great great aunt-uncle",
                          "Great great niece-nephew","First cousin"]:
            return (coef >= 0 and coef <= 0.1)

        elif relation in ["Half sibling","Grandparent","Grandchild","Aunt-Uncle",
                           "Niece-Nehphew"]:
            return (coef >= 0.08838 and coef <= 0.17677)

        elif relation in ["Parent","Child","Sibling","Sibling-other-tissue",
                          "Dizygotic twin","Consanguineous parent"]:
            return (coef >= 0.17677 and coef <= 0.35355)

        elif relation in ["Proband","Proband-other-tissue","Monozygotic twin"]:
            return (coef >= 0.45 and coef <= 0.5)

        else:
            raise Exception("relation {0} not found in check !!".format(relation))    
        
    def update_db(self,statement):
        """
        Update database with the statement 
        """

        print("using '{}'".format(statement))
        return

        db = get_connection()
        try:
            cur = db.cursor()
            cur.execute(statement)
            db.commit()
        finally:
            if db.open:
                db.close()

    def is_repcon_fam_mem(self,sample_name):
        """
        Return true/false depending on whether the given samples is a 
        representative control family member 
        """
        
        db = get_connection()
        try:
            cmd = ("""SELECT RepConFamMem FROM SampleT """
                   """WHERE CHGVID = '{0}' """.format(sample_name))
            cur = db.cursor()
            cur.execute(cmd)
            row = cur.fetchone()
            if row:
                if row[0] == "1":
                    return True
                else:
                    return False
            else: ## Should not happen
                raise Exception("CHGVID not present in SampleT table !!")
        finally:
            if db.open:
                db.close()

    def get_proband(self,family_name):
        """
        """
        
        db = get_connection()
        try:
            cmd = (""" SELECT CHGVID FROM SampleT WHERE FamilyRelationProband = 'Proband' """
                   """ AND FamilyID = '{0}'""".format(family_name))
            cur = db.cursor()
            cur.execute(cmd)
            row = cur.fetchone()
            if row:
                return row[0][0]
            else:
                return None
                print Exception("Proband not found for FamilyID : {0}".format(family_name))
        finally:
            if db.open:
                db.close()
                
        
    def update_within_family(self):
        """ Check and update within family relatedness table
        """

        prev_family_id = None


        statement = ("""INSERT INTO within_fam_relatedness """
                     """(FamilyID, sample_name1, pseudo_prepid1 """
                     """,sample_name2, pseudo_prepid2, coef) """
                     """ VALUES ('{0}','{1}',{2},'{3}',{4},{5}) ON DUPLICATE KEY """
                     """ UPDATE coef = {5} """)
        

        if self.parse_kin_within(self.king_wtn):
            for family_id,sample_name,pseudo,related_sample_name,related_pseudo,coef in self.parse_kin_within(self.king_wtn):
                proband = self.get_proband(family_id)
                ## Try making the code below more elegant , it's a mess now ! 
                if sample_name != proband and related_sample_name != proband:
                    if self.is_repcon_fam_mem(sample_name) and self.is_repcon_fam_mem(related_sample_name):
                        if coef >= 0.1:
                            self.update_db( statement.format( family_id, sample_name, pseudo,related_sample_name, related_pseudo,coef ) )

                elif sample_name == proband:
                    relation = self.get_relation_to_proband(related_sample_name)
                    if not check_relation(relation,coef):
                        self.update_db( statement.format( family_id,sample_name, pseudo,related_sample_name, related_pseudo,coef ) )
                elif proband != None:
                    relation = self.get_relation_to_proband(sample_name)
                    if not check_relation(relation,coef):
                        self.update_db(statement.format(family_id,sample_name, pseudo,related_sample_name, related_pseudo,coef) )
    
                    
                
    def parse_kin_within(self,f):
        """ Parse within family king output 
        """
        
        i = 0
        with open(f,'r') as IN:
            for line in IN:
                line=line.strip('\n')
                contents = line.split('\t')
                if contents[0] == 'FID':
                    continue
                else:
                    i+=1
                    coef = float(contents[-2])
                    family_id = contents[0]
                    sample_name = contents[1].split('_')[0]
                    pseudo = contents[1].split('_')[1]
                    related_sample_name = contents[2].split('_')[0]
                    related_pseudo = contents[2].split('_')[1]
                    yield [family_id,sample_name,pseudo,related_sample_name,related_pseudo,coef]
                

                    
    def parse_kin_btw(self,f):
        """ Parse between family king output 
        """

        i = 0
        with open(f,'r') as IN:
            for line in IN:
                line=line.strip('\n')
                contents=line.split('\t')
                #print contents
                if contents[0] == 'FID1': # header
                    continue
                else:
                    i+=1
                    coef = float(contents[-1])
                    if coef > self.btw_fam_threshold: # sample is above cryptic relatedness threshold
                        sample_name = contents[1].split('_')[0]
                        pseudo = contents[1].split('_')[1]
                        if contents[1] == contents[0]: # if family name is the same as the sample name,
                            ## use the parsed output from above, else use the actual value.
                            family_id = sample_name
                        else:
                            family_id = contents[0]
                            ## Note : The above shenanigan has to be done because sample names have to be unique in the ped file
                            ## defining the sample_name/chgvid in combination with the pseudo_prepid will ensure this
                            ## We could have also just have used the pseudo_prepids (maybe it might be cleaner to code ?)
                            ## The sample names are defined as '{sample_name}_{pseudo_prepid}' in the masterped file

                            ## The family name is defined to be the same as the sample_name if it is not available in SampleT 
                            
                        related_sample_name = contents[3].split('_')[0]
                        related_pseudo = contents[3].split('_')[1]
                        if contents[2] == contents[3]: ## same logic as above
                            related_family_id = related_sample_name
                        else:
                            related_family_id = contents[2]

                        yield [sample_name,pseudo,family_id,related_sample_name,related_pseudo,related_family_id,coef]
                       
        
class RunKing(SGEJobTask):
    """ Run king
    """

    n_cpu = 1
    parallel_env = "threaded"
    shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    masterped = luigi.Parameter(default="/nfs/seqscratch09/rp2801/relatedness_check/masterped.ped")
    output_directory = luigi.Parameter()
    markers = luigi.Parameter(default="{0}/filtered.37MB.master.training.map".format(os.getcwd()))
    genotyping_dir  = luigi.Parameter()
    batch_size = luigi.IntParameter()
        
    def __init__(self,*args,**kwargs):
        super(RunKing,self).__init__(*args,**kwargs)
        self.log_file = os.path.join(self.output_directory,"relatedness_check.log")
        self.temp_bed = os.path.join(self.output_directory,"temp_masterped")
        self.output_prefix = os.path.join(self.output_directory,"relatedness")
        self.king_btw           = self.output_prefix+".kin0"
        self.king_wtn           = self.output_prefix+".kin"
        self.verification_file = os.path.join(self.output_directory,"king.success")
        self.masterped_stem = self.masterped[0:-4]
        self.run_locally=True
        
    def work(self):

        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)

        plink_cmd = (" /nfs/seqscratch11/rp2801/bin_seqscratch/plink --make-bed --file {0} --out {1} &> {2}""".format(
          self.masterped_stem, self.temp_bed, self.log_file ) )
        
        print("\n\nrunning '{}'\n\n".format(plink_cmd))
        run_shellcmd(plink_cmd)
        
        king_cmd = (" /nfs/goldstein/software/sh/king -b {0}.bed --kinship --related --degree 3 --prefix {1} &> {2}""".format(
          self.temp_bed, self.output_prefix, self.log_file ) )

        print("\n\nrunning '{}'\n\n".format(king_cmd))
        run_shellcmd(king_cmd)
        
        verification_cmd = """ touch {0}""".format(self.verification_file)

        run_shellcmd(verification_cmd)
            
    def requires(self):
        ###### this is where 'clone' is usefull - i.e. don't need to re-provide the mutliply specified output_directory=luigi.Parameter when used in both this and dep?!?
        return self.clone(AppendMasterPed)
        # return AppendMasterPed( masterped=self.masterped, output_directory=self.output_directory, markers=self.markers, run_locally=True, genotyping_dir=self.genotyping_dir )

    def output(self):
        return luigi.LocalTarget(self.verification_file)


class AppendMasterPed(SGEJobTask):

    # n_cpu = 1
    # parallel_env = "threaded"
    # shared_tmp_dir = "/nfs/seqscratch09/tmp/luigi_test"

    masterped = luigi.Parameter()
    output_directory = luigi.Parameter()
    markers = luigi.Parameter(default="{0}/filtered.37MB.master.training.map".format(os.getcwd()))
    genotyping_dir  = luigi.Parameter()
    batch_size = luigi.IntParameter()
    
    def __init__(self,*args,**kwargs):

        super(AppendMasterPed,self).__init__(*args,**kwargs)

        self.samples_to_append = get_samples_append()[0:self.batch_size]
        self.log = os.path.join(self.output_directory,"append_masterped.log")
        self.cmd = "cat {0} >> {1} 2>> {2}"
        self.run_locally=True
        
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)
        
    def work(self):
        
        for sample_name, pseudo_prepid, sample_type in self.samples_to_append:

            pedfile = os.path.join(self.genotyping_dir, "{0}.{1}.{2}/{0}.{1}.ped".format( sample_name, pseudo_prepid, sample_type ) )

            run_shellcmd( self.cmd.format( pedfile, self.masterped,self.log ) )

            update_ped_status( pseudo_prepid, "append_ped" )
                    
    def requires(self):
        for chgvid, pseudo, seq in self.samples_to_append:
            # yield self.clone(CreatePed) 
            yield CreatePed( pseudo_prepid=pseudo, sample_name=chgvid, align_loc='irrelevant', sample_type=seq, 
                             output_directory=self.output_directory, markers=self.markers, run_locally=True )
        
    def output(self):
        for chgvid,pseudo,seq in self.samples_to_append:
            yield SQLTarget(pseudo,"append_ped")

  
