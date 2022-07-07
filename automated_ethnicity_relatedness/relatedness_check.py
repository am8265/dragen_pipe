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

################## this entire things seems pointless. running ethnicity as part of pipeline and just pre-check the combined ped for whether something 
################## has been integrated and concatenate any/all missing. the questions are really more about consistency checks, counts for categories,
################## ignoring things for which the value doesn't change e.g. some form of md5sum?!?

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
    masterped       = luigi.Parameter(default="/nfs/fastq_temp2/dsth/relatedness/combined.ped")
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
        # self.btw_fam_threshold  = 0.1
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
        print("\n\n> between\n\n")
        statement = ("INSERT INTO between_fam_relatedness (sample_name1, pseudo_prepid1, FamilyID1,sample_name2, pseudo_prepid2, FamilyID2, coef) \
                     VALUES ('{0}',{1},'{2}','{3}',{4},'{5}',{6}) ON DUPLICATE KEY UPDATE coef = {6} ")
        if self.parse_kin_btw(self.king_btw): ## The file was not empty , an empty file means no relatedness was inferred in

            # why bother naming them?!?
            for sample_name,pseudo,family_name,related_sample_name,related_pseudo,related_family_name,coef in self.parse_kin_btw(self.king_btw):

                self.update_db( statement.format( sample_name,pseudo,family_name,related_sample_name,related_pseudo,related_family_name,coef ) )

######### ONLY RUN WES & WGS!?!?

#### ped format : #Family  Subject  Father  Mother  Sex Phenotype Marker1 M2      Markr3  

####### OLD VERSION SEEMS TO IGNORE FAMILYID 

################ presumably repconfam is the 'unique' control for a family - i.e. not supposed to be closely related to one another e.g. parents?!?
################ presumably repconfam is the 'unique' control for a family - i.e. not supposed to be closely related to one another e.g. parents?!?
################ presumably repconfam is the 'unique' control for a family - i.e. not supposed to be closely related to one another e.g. parents?!?

# from form : Representative 'Control' Family Member : The Project Lead must select one sample per family to be used as a control in other projects; 
# this is not related to affectation status. Control = '1' and all other relatives = '0'. 
# For TRIOS only: label both parents as 1 (not related to each other and thus can both be used as controls).

####### goes through 'within' family (.kin) output : FID     ID1     ID2     N_SNP   Z0      Phi     HetHet  IBS0    Kinship Error                                          
####### THE WITHIN FAMILY FORMAT .kin GIVES ALL PAIR-WISE RELATIONSHIPS WITHIN A FAMILYID i.e. single-member->0,2->1,3->3...

####### checks samples that 'should' be identical i.e. 'chgvid' are within thresholds for 'Proband-other_tissue'... - else is screwed on 'identity' check - i.e. kinship wrt., min/max message
####### else checks if both are 'repconfam' then issues 'repconfam' error if kinship is > 0.1?!? - i.e. they are NEVER supposed to be related!?!
####### else it's 'assumed' same family when one or other is proband - then compares the relationship - via pulling from first and then second individual in sequence
#######     to it's thresholds and gives problem 'relationship' kinship (thresholds)
####### else complains the relationship is unknown

####### goes through 'between' family members .kin0 file - i.e. presumably all pairwise relations between non-family members - i.e. famid1 & famid2 should differ?!?
####### if kingship is above 0.1 (cryptic) then get cryptic-check error with the sample/kinship/threshold - seems NOT to be list i.e. just over-writes errors as opposed to giving count!?!

####### overall has 3 fields RelatednessCheck=Pass/Fail/Not_Checked, FamilyRelatedness=appears-to-be-concatenation-of-relatives, RelatednessError=concatenation-of-problems?!?
####### this is just the expectedrelatedness_06222014.txt info - that ONLY uses cut -f1-3 : 'relationship','min','max'...
####### this is just the expectedrelatedness_06222014.txt info - that ONLY uses cut -f1-3 : 'relationship','min','max'...
####### this is just the expectedrelatedness_06222014.txt info - that ONLY uses cut -f1-3 : 'relationship','min','max'...

# mysql -udh2880 -p1qaz@WSX -h10.73.50.38 sequenceDB -P53306 -e 'select RelatednessCheck,RelatednessError,FamilyRelatedness from seqdbClone where RelatednessError != "None" or FamilyRelatedness != "None"'
# Pass    None    38254.palmfree01.exome:0.5000;
# Pass    None    38585.palmfree06.exome:0.4996;
# Pass    None    38255.palmfree09.exome:0.4997;
# Pass    None    38256.palmfree10.exome:0.5000;
# Pass    None    38350.palmfree12.exome:0.5000;
# Fail    Failed Cryptic Check 31960.epifam01301bip4.exome (0.1206>0.100);        15.irishepiep01303b1.exome:0.0949;
# Fail    Failed Cryptic Check 31960.epifam01301bip4.exome (0.2285>0.100);        14.irishepiep01297b2.exome:0.0949;
#...
# Pass    None    57.isnd26080ac3.exome:0.2532;56.isnd26079ac2.exome:0.2471;
# Pass    None    57.isnd26080ac3.exome:-0.0230;55.isnd25793ac1.exome:0.2471;
# Pass    None    56.isnd26079ac2.exome:-0.0230;55.isnd25793ac1.exome:0.2532;
# mysql -udh2880 -p1qaz@WSX -h10.73.50.38 sequenceDB -P53306 -e 'select RelatednessCheck,RelatednessError,FamilyRelatedness from seqdbClone where RelatednessError != "None" or FamilyRelatedness != "None" and RelatednessCheck != "Pass"' | grep -v 'Failed Cryptic Check' | grep -v 'Not checked'
# Fail    Failed 598.pngn026.exome Proband 0.2815 not (0.4500,0.5000);    598.pngn026.exome:0.2815;
# Fail    Failed 623.pngn0044a1.exome Proband 0.2613 not (0.4500,0.5000); 623.pngn0044a1.exome:0.2613;
# Fail    Failed 622.pngn0045a2.exome Proband 0.2613 not (0.4500,0.5000); 622.pngn0045a2.exome:0.2613;
# Fail    Failed 631.pngn0034.exome Proband 0.2505 not (0.4500,0.5000);   631.pngn0034.exome:0.2505;
# Fail    Failed 630.pngn0035.exome Proband 0.2505 not (0.4500,0.5000);   630.pngn0035.exome:0.2505;
# Fail    Failed 27184.KCTN059e5.exome Sibling -0.0211 not (0.1768,0.3536);Failed 22328.kctn012e3.exome Sibling -0.0397 not (0.1768,0.3536);      27213.KCTN055e1.exome:-0.
# Fail    Failed 27213.KCTN055e1.exome Aunt-Uncle -0.0148 not (0.0884,0.1768);Failed 27194.KCTN057e3.exome Parent 0.0026 not (0.1768,0.3536);Failed 27189.KCTN058e4.exome P
# Fail    Failed 27195.KCTN049c1.exome Proband -0.0085 not (0.4500,0.5000);Failed 27190.KCTN050c2.exome Sibling -0.0126 not (0.1768,0.3536);Failed 22334.kctn003c1.exome Pa
# Fail    Failed 27195.KCTN049c1.exome Sibling 0.0083 not (0.1768,0.3536);Failed 22324.kctn005c3.exome Sibling -0.0102 not (0.1768,0.3536)

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

        # print("using '{}'".format(statement))
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
        print("\n\n> within\n\n")
        prev_family_id = None
        statement = ("INSERT INTO within_fam_relatedness FamilyID, sample_name1, pseudo_prepid1,sample_name2, pseudo_prepid2, coef) \
                     VALUES ('{0}','{1}',{2},'{3}',{4},{5}) ON DUPLICATE KEY UPDATE coef = {5} ")
        if self.parse_kin_within(self.king_wtn):
            for family_id,sample_name,pseudo,related_sample_name,related_pseudo,coef in self.parse_kin_within(self.king_wtn):

                ##### this surely can't be right?!?
                ##### this surely can't be right?!?
                ##### this surely can't be right?!?
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
        print(" > checking file '{}'".format(f))
        i = 0
        with open(f,'r') as IN:
            ### not actually using Phi - i.e. set to 0 - so all will be in error unless not related...?!?
            for line in IN:
                # FID     ID1     ID2     N_SNP   Z0      Phi     HetHet  IBS0    Kinship Error
                line=line.strip('\n')
                contents = line.split('\t')
                if contents[0] == 'FID':
                    continue
                else:
                    i+=1

                    coef = float(contents[-2])
                    print("within family relation '{}'".format(coef))

                    family_id = contents[0]
                    sample_name = contents[1].split('_')[0]
                    pseudo = contents[1].split('_')[1]

                    related_sample_name = contents[2].split('_')[0]
                    related_pseudo = contents[2].split('_')[1]

                    yield [family_id,sample_name,pseudo,related_sample_name,related_pseudo,coef]
                

                    
    def parse_kin_btw(self,f):
        print(" > checking file '{}'".format(f))
        i = 0
        with open(f,'r') as IN:
            for line in IN:
                line=line.strip('\n')
                # FID1    ID1     FID2    ID2     N_SNP   HetHet  IBS0    Kinship                                                        
                contents=line.split('\t')
                if contents[0] == 'FID1': # header
                    continue
                else:
                    i+=1
                    coef = float(contents[-1])
                    if coef > 0.1: 
                        print("they're related '{}'".format(coef))

                        sample_name = contents[1].split('_')[0]
                        pseudo = contents[1].split('_')[1]

                        ##### wtf?!? : why not just treat everything systematically instead?!?
                        # if family name is the same as the sample name, use the parsed output from above, else use the actual value.
                        if contents[1] == contents[0]: 
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

  
