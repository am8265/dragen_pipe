#!/usr/bin/env python2.7
import os,sys,string,subprocess,pdb,re
# import globals

#150127 slg added DBhost command line arg.

# DBuser     = globals.DBUSER
# DBpassword = globals.DBPW
# DB         = globals.DBNAME


print 'STARTED RUNNING RELATEDNESS CHECK'

#if len(sys.argv) != 6:
    # print sys.argv[0], ": wrong number of command line args.  Exiting." 
    # sys.exit(1)
scratch = 'argh';
##### king authors give guidelines for these?!?
thresholdFile = 'expectedrelatedness_06222014.txt'
cryptic_threshold = 1.0
# DBhost = sys.argv[4]                    #provide db host ip
# testMode = int(sys.argv[5])             #This mode uses the Bioifo database and not the main SeqDB

RUNdir = scratch + 'NEWTEST'              #This is the directory where KING runs
plinkOUT = scratch + 'NEWTEST/newbedout'  #Name of Plink PED/MAP binary file that KING uses
plink_software = '/nfs/goldstein/goldsteinlab/software/PLINK_107/plink'  
king_software = '/nfs/goldstein/goldsteinlab/software/king_relatedness/king'
master_ped = '/nfs/informatics/data/pipeline/dragen_pipe/RELATEDNESS_QC/test'
                                        
'''
All threshold based on kinship values from KING
KINGs manual - http://people.virginia.edu/~wc9c/KING/manual.html
KING ouputs two files - king.kin and king.kin0. king.king has the relatedness information wihin family for all samaples and king.kin0 has kinship for unrelated samples
'''

# line = plink_software + " --file " + master_ped + " --make-bed --out " + plinkOUT + " --noweb"
# line = king_software + " -b newbedout.bed --kinship --related --degree 3"

# Reading in threshold file that defines kinship coefficients for a given relationship 
thresholds_min = {}
thresholds_max = {}

inFile = open(thresholdFile,'r')
for line in inFile:
    temp = line.split("\t")
    relationship = temp[0]
    min = temp[1]
    max = temp[2]
    thresholds_min[relationship]=min
    thresholds_max[relationship]=max

from pprint import pprint as pp

# This function reads in a PED id (3751.van2575.genome) and outputs the prepID
def parsePedID(pedID):
    array = None
    prepID = None
    array = pedID.split(".")
    prepID = array[0]

    if prepID.isdigit():
        return prepID
    else:
        return "-9"

# This function reads in a PED id (3751.van2575.genome) and outputs the seqType
def parseSeqType(pedID):
    array = None
    seqType = None
    array = pedID.split(".")
    seqType = array[2]
    
    if seqType:
        return seqType
    else:
        return "-9"
    

def getPrepID(): # Get all unique PrepIDs
 line = "mysql -u " + DBuser + "  -p" + DBpassword + " -h " + DBhost + " -e \"use " + DB + "; select prepID,chgvID,seqType from seqdbClone\" "
 tempPrepIDs = []
 tempID = subprocess.Popen(line, stdout = subprocess.PIPE, shell=True)

 sqlout = tempID.stdout.read().split("\n")
 sqlout.pop()
 sqlout.pop(0)

 Hash={}

 for s in sqlout:
     
     ###Getting the sequence type to be lower case
     id = None
     seqType = None
     array = s.split("\t");
     seqType = array[2].lower()
     id = array[0]+ "." + array[1] + "." + seqType
     Hash[id]=1

 if len(Hash.keys()) > 0 :    
     return Hash
 else:
     return "-9"


def getFamID():# Return all family IDs from database
    
    line =  "mysql -u " + DBuser + "  -p" + DBpassword + " -h " + DBhost + " -e \"use " + DB + "; select distinct FamilyID from seqdbClone\" "
    tempFamilies = []
    tempF = subprocess.Popen(line, stdout = subprocess.PIPE, shell=True)
  
    sqlout = tempF.stdout.read().split("\n") 
    sqlout.pop()
    sqlout.pop(0)
      
    if len(sqlout) > 0 :
        return sqlout
    else:
        return "-9"

def getFamilyMembers(famID):#Given a family ID, return all related individulas if family has more than 1 member
     line = "mysql -u " + DBuser + "  -p" + DBpassword + " -h " + DBhost + " -e \"use " + DB + "; select prepID,chgvID,seqType,RepConFamMem,FamilyRelationProband from seqdbClone where familyID = '"+famID+"' \" "
     tempRelated = []
     tempR = subprocess.Popen(line, stdout = subprocess.PIPE, shell=True)

     sqlout = tempR.stdout.read().split("\n")
     sqlout.pop()
     sqlout.pop(0)
     
     if len(sqlout) > 1 :
         return sqlout
     else:
         return "-9"
     
     

def getRepConFam(famID):#Given a family ID, return all RepConFamMem with values == 1

  line = "mysql -u " + DBuser + "  -p" + DBpassword + " -h " + DBhost + " -e \"use " + DB + "; select prepID from seqdbClone where RepConFamMem = '1' and familyID = '"+famID+"' \" "

  tempsamples = []
  tempx = subprocess.Popen(line, stdout = subprocess.PIPE, shell=True)
  
  sqlout = tempx.stdout.read().split("\n") # this prints something like this ['prepID', '1385', '3772', '3802', '']
  #So I will remove the first and last items from that array
  sqlout.pop()
  sqlout.pop(0)
  
  if len(sqlout) > 0 :
      return sqlout
  else:
      return "-9"



def getProband(famID):#Given a family ID, return all Probands

      line = "mysql -u " + DBuser + "  -p" + DBpassword + " -h " + DBhost + " -e \"use " + DB + "; select prepID from seqdbClone where RepConFamMem = '1' and familyID = '"+famID+"' \" "

      tempsamples = []
      tempx = subprocess.Popen(line, stdout = subprocess.PIPE, shell=True)
      sqlout = tempx.stdout.read().split("\n") # this prints something like this ['prepID', '1385', '3772', '3802', '']
      
      #So I will remove the first and last items from that array
      sqlout.pop()
      sqlout.pop(0)

      if len(sqlout) > 0 :
          return sqlout
      else:
          return "-9"
      
'''
This function queries seqDB and gets the relationship to proband
'''
def getProbandRelationship(sampleid,seqtype,prepID):# return relationship to proband info
    
    line = "mysql -u " + DBuser + "  -p" + DBpassword + " -h " + DBhost + " -e \"use " + DB + "; select FamilyRelationProband from seqdbClone where prepID = "+prepID+" and chgvid =\"" + sampleid + "\" and seqtype = \""+ seqtype +"\";\'"
    
    tempsamples = []
    tempx = subprocess.Popen(line, stdout = subprocess.PIPE, shell=True)
    sqlout = tempx.stdout.read().split("\n")

    if len(sqlout) > 1 :
        return sqlout[1].strip("\n")
    else:
        return "-9"



'''
This function checks to see if two PED ids of the form
PREPID.chgvID.seqtype have the same CHGVID
'''
    
def checkIdentity(sample1, sample2):
    identity = None
    array = None
    array = sample1.split(".")
    chgvID_sample1 = array[1]
    
    array = None
    array = sample2.split(".")
    chgvID_sample2 = array[1]

    if chgvID_sample1 == chgvID_sample2:
        identity = 1
    else:
        identity = 0
        
    return identity

'''
This function appends a value to an existing key value
'''
        
def updateHashMessage(key,message,hash):
    if key in hash.keys():
        str1 = hash[key]
        new_message = message + str1
    else:
        new_message = message
        
    return new_message
    



## Get a list of all samples in Sequence

### READING IN PED FILE
def ReadPED(): # Create a Hash of IDs with PED file IDs
    hash={}
    pedfile = open(scratch+"test.ped")
    for line in pedfile:
        temp = line.split(" ")
        sample = temp[1]
        hash[sample]=1

    return hash


## First populate hashes
prepIDHash = {}
prepIDHash = getPrepID()
PEDHash = ReadPED()

#QC specific Hash
QC_Fam_Relatedness = {} # get kinship for all family members 
QC_Fail_message = {}    # acutal message for each sample, append to QC_FAM_Relatedness
QC_Check = {}           # Pass/Fail/Not Checked     

# Other Hashes which will be initialized below
ProbandHash = {}                      # contains probands, as specified by SeqDB
RelationshipProbandHash = {}          # specifies the relationship to proband
RepConFamHash = {}                    # determines if samples is 0/1

#Get string containing all distinct family members with greater than 1 individual
family = getFamID()
for f in family:
    familyList = getFamilyMembers(f)
    if familyList != '-9':

        for famMember in familyList:
            array = famMember.split("\t")
            prepID = array[0]
            chgvID = array[1]
            seqType = array[2].lower()
            repConFamMember=array[3]
            relationProband=array[4]

            pedID = prepID + '.' + chgvID + '.' + seqType
            
            ### Populating other hashes
            RelationshipProbandHash[pedID] = relationProband

            if repConFamMember == '1':
                RepConFamHash[pedID] = 1

            if relationProband == 'Proband':
                ProbandHash[pedID] = 1


## Opening king file
os.chdir(RUNdir)
kinfile = open("king.kin").readlines() #these are all individuals in a family

for line in kinfile[1:]:
    sample1=None
    sample2=None
    kinship=None
    min=None
    max=None
    relationship=None
    array=[]
    str1=None
    str2=None
    
    line = line.strip("\n")
    temp = line.split("\t")
    sample1 = temp[1]
    sample2 = temp[2]
    kinship = temp[-2]
    kinship_numeric = float(temp[-2])

    # Populating the family relatedness hash
    list1 = [sample1,sample2]
    list2 = [sample2,sample1]
    
    for i in range (len(list1)):
        s1 = list1[i]
        s2 = list2[i]

        message = s2 + ':' + kinship + ';'

        new_message = updateHashMessage(s1,message,QC_Fam_Relatedness)
        QC_Fam_Relatedness[s1] = new_message


    ### Check for identity
    if checkIdentity(sample1,sample2) == 1:
        max = float(thresholds_max['Proband-other tissue'])
        min = float(thresholds_min['Proband-other tissue'])

        min_string = "%.4f" % min
        max_string = "%.4f" % max

        ###### max with identity check?!?
        if kinship_numeric < min or kinship_numeric > max:

            list1 = [sample1,sample2]
            list2 = [sample2,sample1]
            
            for i in range (len(list1)):
                s1 = list1[i]
                s2 = list2[i]
                message = 'Failed Identity' + s2 + ' ' + kinship + ' not (' + min_string + ',' + max_string + ');'
                new_message = updateHashMessage(s1,message,QC_Fail_message)
                QC_Fail_message[s1] = new_message
                    
    ## If both samples are designated as RepConFamily Member ==1 then have to check threshold        
    elif (sample1 in RepConFamHash.keys() and sample2 in RepConFamHash.keys()):  
        if kinship_numeric > cryptic_threshold:
            cryptic_threshold_string = "%.3f" % cryptic_threshold

            list1 = [sample1,sample2]
            list2 = [sample2,sample1]
            for i in range (len(list1)):
                s1 = list1[i]
                s2 = list2[i]
                message = 'Failed RepConFamMember ' + s2 + ' (' + kinship + '>' + cryptic_threshold_string + ');'                        
                new_message = updateHashMessage(s1,message,QC_Fail_message)                    
                QC_Fail_message[s1] = new_message

    #### Checking relationship between samples in the same family         
    elif (sample1 in ProbandHash.keys() or sample2 in ProbandHash.keys()):  #### Checking relationship between samples
        if sample1 in ProbandHash.keys():
            relationship = RelationshipProbandHash[sample2]
        elif sample2 in ProbandHash.keys():
            relationship = RelationshipProbandHash[sample1]

        if relationship in thresholds_max and relationship in thresholds_min:
            max = float(thresholds_max[relationship])
            min = float(thresholds_min[relationship])
        
            if kinship_numeric < min or kinship_numeric > max:
                min_string = "%.4f" % min
                max_string = "%.4f" % max

                list1 = [sample1,sample2]
                list2 = [sample2,sample1]
                for i in range (len(list1)):
                    s1 = list1[i]
                    s2 = list2[i]
                    message = 'Failed ' + s2 + " " + relationship + " "+ kinship + ' not (' + min_string + ',' + max_string + ');'
                    new_message = updateHashMessage(s1,message,QC_Fail_message)
                    QC_Fail_message[s1] = new_message
                
        else:
            print 'Min, Max  kinship coefficient not defined for relationship', relationship, sample1, sample2

### NOW opening the Cryptic relatedness fails
os.chdir(RUNdir)
kinfile = open("king.kin0").readlines() #these are all individuals in a family

for line in kinfile[1:]:
    sample1=None
    sample2=None
    kinship=None
    kinship_numeric=None

    line = line.strip("\n")
    temp = line.split("\t")
    sample1 = temp[1]
    sample2 = temp[3]
    kinship = temp[7]
    kinship_numeric = float(temp[7])

    if kinship_numeric > cryptic_threshold:
        cryptic_threshold_string = "%.3f" % cryptic_threshold

        list1 = [sample1,sample2]
        list2 = [sample2,sample1]
        for i in range (len(list1)):
            s1 = list1[i]
            s2 = list2[i]
            message = 'Failed Cryptic Check ' + s2 + ' (' + kinship + '>' + cryptic_threshold_string + ');'
            new_message = updateHashMessage(s1,message,QC_Fail_message)                
            QC_Fail_message[s1] = new_message
            # print 'Failed Cryptic Check ', message_sample
    
 ### END of cryptic relatedness file                       


### Final loop
for sample in prepIDHash.keys(): # this is a list of all samples in the database with prepID.chgvID.seqType
    sample_seqType = parseSeqType(sample)
    sample_seqType = sample_seqType.lower()
  #If sample is custom capture or RNA-seq enter relatedness check as N/A
    if  re.search('custom_capture',sample_seqType) or re.search('rnaseq',sample_seqType):
        QC_Check[sample]='NA'
        QC_Fam_Relatedness[sample] = 'NA'
        QC_Fail_message[sample] = 'NA'
      
    elif sample not in PEDHash.keys():  # this is a list of samples in the PED file
        QC_Check[sample]='Not checked'
        QC_Fam_Relatedness[sample] = 'NA'
        QC_Fail_message[sample] = 'NA'       
    else: #sample in PED hash
        if sample in QC_Fail_message.keys():
            QC_Check[sample]='Fail'
        else: # Pass QC
            QC_Check[sample]='Pass'
            QC_Fail_message[sample] = 'None';
          
        if sample not in QC_Fam_Relatedness.keys():
            QC_Fam_Relatedness[sample] = 'None';
              
## Update the database
            
for sample in prepIDHash.keys():
    prepID          = parsePedID(sample)
    check           = QC_Check[sample]
    fail_message    = QC_Fail_message[sample]
    fam_relatedness = QC_Fam_Relatedness[sample]
    
    # slg set RelatednessNotes (deprecated) to NULL 14/10/30 to delete old inconsistent values
    line = "mysql -u " + DBuser + "  -p" + DBpassword + " -h " + DBhost + " -e \"use " + DB +"; "
    line = line + " update seqdbClone set RelatednessCheck = '"+check+"',RelatednessError = '"+fail_message+"',FamilyRelatedness = '"+fam_relatedness+"', RelatednessNotes=NULL "
    line = line + " where prepID = '"+prepID+"' \" "
    #line = "mysql -u " + DBuser + "  -p" + DBpassword + " -h " + DBhost + " -e \"use " + DB + "; update seqdbClone set RelatednessCheck = '"+check+"',RelatednessError = '"+fail_message+"',FamilyRelatedness = '"+fam_relatedness+"' where prepID = '"+prepID+"' \" "
    os.system(line)


print 'FINISHED RUNNING RELATEDNESS CHECK'
