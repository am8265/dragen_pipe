#!/usr/bin/env python
import sys
from optparse import OptionParser 
import subprocess
import os 
import regionstobed 


def get_filename_from_fullpath(full_path):
    ''' Get just the filename from a full path file
    '''
    
    return full_path.split('/')[-1]


def run_cmd(cmd_to_run):
    ''' Run a command
    '''
    os.system(cmd_to_run) ## Replace this with subprocess 

    
def run_markduplicates(options):
    ''' Run Markduplicates
    '''
    ## Unpack arguments
    java_location = options.java_location
    picard_location = options.picard_location
    bam_file = options.bam_file 
    output_file = options.output_dir + options.sample +".markduplicates.txt"
    temp_output_bam = options.scratch_dir + 'temp.bam'

    ## Define the picard command to run 
    markduplicates_cmd = "%s -jar %s MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s"%(java_location,picard_location,bam_file,temp_output_bam,output_file)
    ## Execute the command
    run_cmd(markduplicates_cmd)

    ## Remove the temp bam file
    os.system('rm %s'%(temp_output_bam))

    
def create_targetfile(options):
    ''' Use this function to create a target interval file list
    compatible with GATK and Picard Tools
    '''
    ## Unpack arguments for BedToIntervalList
    java_location = options.java_location
    picard_location = options.picard_location
    bed_file = options.bed_file 
    seqdict_file = options.seqdict_file

    
    output_name = get_filename_from_fullpath(bed_file)
    output_file = options.output_dir + output_name + ".list"
  
    bedtointerval_cmd = "%s -jar %s BedToIntervalList I=%s SD=%s OUTPUT=%s"%(java_location,picard_location,bed_file,seqdict_file,output_file)
    run_cmd(bedtointerval_cmd)

    
def run_wgsmetrics(options):
    ''' Run Picard's CollectWgsMetrics Function for
    Genome sequencing samples
    '''

    ## Unpack arguments
    java_location = options.java_location
    picard_location = options.picard_location
    reference_file = options.reference_file
    bam_file = options.bam_file 
    output_file = options.output_dir + options.sample +".metrics.txt"
    
    
    ## Define the CollectWgsMetrics Command
    if(options.wgsinterval == True):
        if(options.create_targetfile == True):
            bed_file = options.bed_file
            output_name = get_filename_from_fullpath(bed_file)
            target_file = options.output_dir + output_name + ".list"
        else:    
            target_file = options.target_file
            
        wgs_cmd = "%s -jar %s CollectWgsMetrics R=%s O=%s I=%s INTERVALS=%s MQ=20 Q=10"%(java_location,picard_location,reference_file,output_file,bam_file,target_file)
        print wgs_cmd
    else:
        wgs_cmd = "%s -jar %s CollectWgsMetrics R=%s O=%s I=%s MQ=20 Q=10"%(java_location,picard_location,reference_file,output_file,bam_file)

    ## Run the command
    run_cmd(wgs_cmd)

    
def run_hsmetrics(options):
    ''' RUN HsMetrics 
    '''
    ## Unpack arguments for hsmetrics
    java_location = options.java_location
    picard_location = options.picard_location
    bait_file = options.bait_file
    target_file = options.target_file
    bam_file = options.bam_file
    sample = options.sample
    output_file = options.scratch + options.sample + ".hsmetrics.txt" 

    if(options.create_targetfile == True):
        bed_file = options.bed_file 
        output_name = get_filename_from_fullpath(bed_file)
        target_file = options.output_dir + output_name + ".list"
    else:    
        target_file = options.target_file

    ## Define the command to run
    hsmetrics_cmd = "%s -jar %s CalculateHsMetrics BI=%s TI=%s METRIC_ACCUMULATION_LEVEL=ALL_READS I=%s O=%s MQ=20 Q=10"%(java_location,picard_location,bait_file,target_file,bam_file,output_file)
    
    ## Run the command
    run_cmd(hsmetrics_cmd)
    

def process_hsmetrics(options):
    ''' Process the HsMetrics Output
    '''

    i=0
    metrics_file = options.scratch + options.sample + ".hsmetrics.txt"
    metrics = []
    headers = []
    with open(metrics_file) as IN:
        for line in IN:
            if(line[0]!='#'):
                i+=1
                if(i==2): # First line after the ## header lines is the headers for our metrics
                    headers=line.strip('\n').split('\t')
                if(i==3): # The metrics line is the 2nd line after the ## header lines  
                    metrics=line.strip('\n').split('\t')

    hsmetrics = dict(zip(headers,metrics))

    ##Output the required columns 


def process_wgsmetrics(options):
    ''' Process the wgsMetrics Output
    '''

    i=0
    metrics_file = options.output_dir + options.sample + ".metrics.txt"
    metrics = []
    headers = []
    with open(metrics_file) as IN:
        for line in IN:
            if(line[0]!='#'):
                i+=1
                if(i==2): # First line after the ## header lines is the headers for our metrics
                    headers=line.strip('\n').split('\t')
                if(i==3): # The metrics line is the 2nd line after the ## header lines  
                    metrics=line.strip('\n').split('\t')
    
    
    wgsmetrics = dict(zip(headers,metrics))
    ## Output required columns (PCT_1X,PCT_5X,PCT_10X,PCT_15X,PCT20X,MEAN_COVERAGE)

    
def create_bed(options):
    ''' Uses the regionstobed function to create a bed file
    naming scheme for the output bed is simply to append a .bed
    to the regions file name, output is in the same directory as regions file
    Regions file is assumed to have the absolute path, returns the full path to 
    the bed file
    '''

    output_bed_file = options_regions_file + ".bed"
    regionstobed.output_bed(options.regions_file,output_bed_file)
    
    return output_bed_file  


def main(options):
    ''' Main Function
    '''
   
    if(options.create_targetfile == True):
        print "\nCreating The Target Intervals File....\n"
        if options.bed_file == "" and len(options.regions_file > 1): ## If the regions file is used
            options.bed_file = create_bed(options)
        create_targetfile(options)
        print "\nCreated TargetFile\n"
        
    if(options.seq_type.upper() == 'GENOME'):
        print "\nRunning wgsmetrics as this is a GENOME sample....\n"
        if options.wgsinterval == True:
            if options.bed_file == "" and len(options.regions_file > 1): ## If the regions file is used
                options.bed_file = create_bed(options) 
                            
        run_wgsmetrics(options)
        process_wgsmetrics(options)
        print "\nFinished Running Metrics\n"
        
    elif(options.seq_type.upper()!= ''):
        print "\nRunning hsmetrics as this is a EXOME/CUSTOM_CAPTURE....\n"
        run_hsmetrics(options)
        process_hsmetrics(options)
        print "\nFinished Running Metrics\n"
        
    
if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("--scratch", type = "string",help="Scratch directory", dest = "scratch", default = ".")
    parser.add_option("--bamfile", type = "string",help="BAM file",dest = "bam_file",default = "")
    parser.add_option("--bed", type = "string",help="BED file", dest = "bed_file",default = "")
    parser.add_option("--target_file", type = "string",help="Target List File(Specify the flag to create it if you dont have one)", dest = "target_file", default = "")
    parser.add_option("--bait_file", type = "string",help="Bait File", dest = "bait_file", default = "")
    parser.add_option("--reference", type = "string",help="Reference Fasta", dest = "reference_file", default = "/nfs/goldsteindata/refDB/HS_Build37/BWA_INDEX_hs37d5/hs37d5.fa")
    parser.add_option("--seq_dictionary", type = "string",help="Sequence Dictionary",dest = "seqdict_file",default="/nfs/goldsteindata/refDB/HS_Build37/BWA_INDEX_hs37d5/hs37d5.dict")
    parser.add_option("--sample", type = "string",help="SAMPLE ID",dest = "sample",default="")
    parser.add_option("--seq_type", type = "string",help="SEQ TYPE : Genome/Exome",dest = "seq_type",default="")
    parser.add_option("--pid", type = "string",help="PREP ID",dest = "prepid",default="")
    parser.add_option("--output",type="string",dest = "output_dir",default = ".")
    parser.add_option("--cluster",type="string",dest = "cluster",default="lsrc")
    parser.add_option("--create_targetfile",action="store_true",help="FLAG-True/False",dest="create_targetfile",default="False")
    parser.add_option("--picard_location",type="string",help="Picard jar full path",dest="picard_location",default="/nfs/goldstein/software/picard-tools-2.2.1/picard.jar")
    parser.add_option("--java_location",type="string",help="Java full path",dest="java_location",default="/nfs/goldstein/software/jdk1.8.0_05/bin/java")
    parser.add_option("--wgsinterval",action="store_true",help="FLAG-True/False",dest="wgsinterval",default="False")
    parser.add_option("--regions_file",action="string",help="ATAV Regions File",dest="regions_file",default="")
    (options,args) = parser.parse_args()
           

    if (options.scratch == "."):
        print "No given scratch directory - using current directory"
    if (options.pythondir == "."):
        print "No given python script directory - using current directory"
    if options.scratch[-1] != "/":
        options.scratch = options.scratch + "/"
    if options.output_dir[-1] != "/":
        options.output_dir = options.output_dir + "/"   
    if (options.bam_file == "" and options.seq_type!=""):
        print "No BAM file - exiting..."
        sys.exit()   
    if (options.bed_file == "" and options.create_targetfile==True and options.regions_file == ""):
        print "No BED file or Regions File provided to create TargetList File - exiting ..."
        sys.exit()
    if (options.target_file == "" and options.create_targetfile==False and options.seq_type.upper()!="GENOME"):
        print "No Target File Provided - exiting..."
        sys.exit()    
    if (options.bait_file == "" and options.seq_type.upper() not in ['GENOME','']):
        print "No Bait File Provided - exiting..."
        sys.exit()
    if (options.sample == "" and options.seq_type.upper()!=''):
        print "No sample id - Exiting..."
        sys.exit()
    if (options.wgsinterval == True and options.bed_file == "" and options.regions_file == ""):
        print "No region specified for wgsmetrics\nTurn of flag if the whole genome metrics are required - Exiting..."
        sys.exit()
    if (options.prepid == ""):
        print "No Prep ID"
    if (options.seq_type == ""):
        print "No SEQ type - No Metrics will be run"
    if (options.picard_location == "/nfs/goldstein/software/picard-tools-2.2.1/picard.jar"):
        print "Using default Picard location : %s"%(options.picard_location)
    if (options.java_location == "/nfs/goldstein/software/jdk1.8.0_05/bin/java"):
        print "Using default Java location : %s"%(options.java_location)
    
      
    main(options)
    
