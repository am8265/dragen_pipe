#!/usr/bin/env Python
import sys
from optparse import OptionParser 
import subprocess
import os 


def get_filename_from_fullpath(full_path):
    ''' Get just the filename from a full path file
    '''
    
    return full_path.split('/')[-1]

def create_targetfile(options):
    ''' Use this function to create a target interval file list
    compatible with GATK and Picard Tools
    '''
    ## Unpack arguments for BedToIntervalList
    picard_location = options.picard_location
    bed_file = options.bed_file 
    seqdict_file = options.seqdict_file

    
    output_name = get_filename_from_fullpath(bed_file)
    output_file = options.output_dir + output_name + ".list"
  
    bedtointerval_cmd = "java -jar %s BedToIntervalList I=%s SD=%s OUTPUT=%s"%(picard_location,bed_file,seqdict_file,output_file)
    os.system(bedtointerval_cmd)

def run_wgsmetrics(options):
    ''' Run GATK's DepthOfCoverage Function for
    Genome sequencing samples
    '''

    ## Unpack arguments for DepthOfCoverage
    picard_location = options.picard_location
    reference_file = options.reference_file
    bam_file = options.bam_file 
    output_file = options.output_dir + options.sample +".metrics.txt"
    java_location = options.java_location 

    depthofcoverage_cmd = "%s -jar %s CollectWgsMetrics -R=%s -O=%s -I=%s MQ=20 Q=10"%(java_location,picard_location,reference_file,output_file,bam_file,target_file)
    os.system(depthofcoverage_cmd)

def run_hsmetrics(options):
    ''' RUN HsMetrics 
    '''
    ## Unpack arguments for hsmetrics
    picard_location = options.picard_dir
    bait_file = options.bait_file
    target_file = options.target_file
    bam_file = options.bam_file
    sample = options.sample
    output_file = options.scratch + options.sample + ".hsmetrics.txt" 

    ## Define the command to run
    hsmetrics_cmd = "java -jar %s CalculateHsMetrics BI=%s TI=%s METRIC_ACCUMULATION_LEVEL=ALL_READS I=%s O=%s"%(picard_dir,baits_file,targets_file,bam_file,output_dir)

    
    ## Run the command 
    
def process_metricsfile():
    ''' Process the HsMetrics Output
    '''

def output_metrics():
    ''' Output the processed file
    '''
    
def main(options):
    ''' Main Function
    '''

    
    if(options.create_targetfile == True):
        create_targetfile(options)


    if(options.seq_type.upper() == 'GENOME'):
        run_depthofcoverage(options)
    
    elif(options.seq_type.upper() == 'EXOME'):
        run_hsmetrics(options)

    #run_hsmetrics()

    #process_metricsfile()

    #output_metrics()

    
if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("--scratch", type = "string",help="Scratch directory", dest = "scratch", default = ".")
    parser.add_option("--bamfile", type = "string",help="BAM file",dest = "bam_file",default = "")
    parser.add_option("--bed", type = "string",help="BED file", dest = "bed_file",default = "")
    parser.add_option("--target", type = "string",help="Target List File(Specify the flag to create it if you dont have one)", dest = "target_file", default = "")
    parser.add_option("--reference", type = "string",help="Reference Fasta", dest = "reference_file", default = "")
    parser.add_option("--seq_dictionary", type = "string",help="Sequence Dictionary",dest = "seqdict_file",default="")
    parser.add_option("--sample", type = "string",help="SAMPLE ID",dest = "sample",default="")
    parser.add_option("--seq_type", type = "string",help="SEQ TYPE : Genome/Exome",dest = "seq_type",default="")
    parser.add_option("--pid", type = "string",help="PREP ID",dest = "prepid",default="")
    parser.add_option("--pythondir", type = "string",help ="Python script directory",dest= "pythondir",default = ".")
    parser.add_option("--output",type="string",dest = "output_dir",default = ".")
    parser.add_option("--cluster",type="string",dest = "cluster",default="lsrc")
    parser.add_option("--create_targetfile",action="store_true",help="FLAG,True/False",dest="create_targetfile",default="False")
    parser.add_option("--picard_location",type="string",help="Picard jar full path",dest="picard_location",default="")
    parser.add_option("--gatk_location",type="string",help="GATK jar full path",dest="gatk_location",default="")
    
    (options,args) = parser.parse_args()

    #if(len(options)<2):
        #print "Specify Arguments ! exiting...\n Use -h/--help for the help menu"
        #sys.exit()

    if (options.scratch == "."):
        print "No given scratch directory - using current directory"
    if (options.pythondir == "."):
        print "No given python script directory - using current directory"
    if options.scratch[-1] != "/":
        options.scratch = options.scratch + "/"
    if options.pythondir[-1] != "/":
        options.pythondir = options.pythondir + "/"
    if options.output_dir[-1] != "/":
        options.output_dir = options.output_dir + "/"   
    if (options.bam_file == "" and options.create_targetfile==False):
        print "No BAM file - exiting..."
        sys.exit()   
    if (options.bed_file == "" and options.create_targetfile==True):
        print "No BED file provided to create TargetList File - exiting ..."
        sys.exit()
    if (options.target_file == "" and options.create_targetfile==False):
        print "No Target File Provided - exiting..."
        sys.exit()    
    if (options.sample == "" and options.create_targetfile==False):
        print "No sample id - Exiting..."
        sys.exit()
    if (options.prepid == ""):
        print "No Prep ID"
    if (options.seq_type == ""):
        print "No SEQ type - No Metrics will be run"
    if (options.picard_location == ""):
        print "Specify Picard Location !\nexiting ..."
        sys.exit()
    if (options.gatk_location == ""):
        print "Specify GATK location !\nexiting ..."
        sys.exit()
    if (options.seqdict_file == "" and create_targetfile == True):
        print "Specify the Sequence Dictionary File To Use!\nexiting ..."
        sys.exit()

    main(options)
    
