#!/usr/bin/env Python
import sys
from optparse import OptionParser 
import subprocess
import os 


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
    metrics_file = options.scratch + options.sample + ".metrics.txt"
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
    
    
def main(options):
    ''' Main Function
    '''
   
    if(options.create_targetfile == True):
        create_targetfile(options)
        print "Created TargetFile\n"
        
    if(options.seq_type.upper() == 'GENOME'):
        run_wgsmetrics(options)
        process_wgsmetrics(options)
    
    else:
        run_hsmetrics(options)
        process_hsmetrics(options)
        
    
if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("--scratch", type = "string",help="Scratch directory", dest = "scratch", default = ".")
    parser.add_option("--bamfile", type = "string",help="BAM file",dest = "bam_file",default = "")
    parser.add_option("--bed", type = "string",help="BED file", dest = "bed_file",default = "")
    parser.add_option("--target_file", type = "string",help="Target List File(Specify the flag to create it if you dont have one)", dest = "target_file", default = "")
    parser.add_option("--bait_file", type = "string",help="Bait File", dest = "bait_file", default = "")
    parser.add_option("--reference", type = "string",help="Reference Fasta", dest = "reference_file", default = "")
    parser.add_option("--seq_dictionary", type = "string",help="Sequence Dictionary",dest = "seqdict_file",default="")
    parser.add_option("--sample", type = "string",help="SAMPLE ID",dest = "sample",default="")
    parser.add_option("--seq_type", type = "string",help="SEQ TYPE : Genome/Exome",dest = "seq_type",default="")
    parser.add_option("--pid", type = "string",help="PREP ID",dest = "prepid",default="")
    parser.add_option("--pythondir", type = "string",help ="Python script directory",dest= "pythondir",default = ".")
    parser.add_option("--output",type="string",dest = "output_dir",default = ".")
    parser.add_option("--cluster",type="string",dest = "cluster",default="lsrc")
    parser.add_option("--create_targetfile",action="store_true",help="FLAG-True/False",dest="create_targetfile",default="False")
    parser.add_option("--picard_location",type="string",help="Picard jar full path",dest="picard_location",default="")
    parser.add_option("--java_location",type="string",help="Java full path",dest="java_location",default="")
    parser.add_option("--wgsinterval",action="store_true",help="FLAG-True/False",dest="wgsinterval",default="False")

    (options,args) = parser.parse_args()
        
        

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
    if (options.bam_file == "" and options.seq_type!=""):
        print "No BAM file - exiting..."
        sys.exit()   
    if (options.bed_file == "" and options.create_targetfile==True):
        print "No BED file provided to create TargetList File - exiting ..."
        sys.exit()
    if (options.target_file == "" and options.create_targetfile==False):
        print "No Target File Provided - exiting..."
        sys.exit()    
    if (options.bait_file == "" and options.seq_type.upper() not in ['GENOME','']):
        print "No Bait File Provided - exiting..."
        sys.exit()
    if (options.sample == "" and options.seq_type.upper()!=''):
        print "No sample id - Exiting..."
        sys.exit()
    if (options.prepid == ""):
        print "No Prep ID"
    if (options.seq_type == ""):
        print "No SEQ type - No Metrics will be run"
    if (options.picard_location == ""):
        print "Specify Picard Location !\nexiting ..."
        sys.exit()
    if (options.java_location == ""):
        print "Specify Java location !\nexiting ..."
        sys.exit()
    if (options.seqdict_file == "" and create_targetfile == True):
        print "Specify the Sequence Dictionary File To Use For Creating the TargetFile!\nexiting ..."
        sys.exit()

    main(options)
    
