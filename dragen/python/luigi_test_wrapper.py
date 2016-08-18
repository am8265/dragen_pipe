#!/nfs/goldstein/software/python2.7.7/bin/python

import glob 
import os
import subprocess
import time 
import sys
import shlex
import argparse

"""
def run_sleep(period):
    proc = subprocess.Popen(['sleep',str(period)])
    return proc

def print_shit():
    proc = subprocess.Popen(['echo','shit'],stdout=subprocess.PIPE)
    return proc

start = time.time()
procs = []
    

for proc in procs:
    out,err = proc.communicate()
    print out

end = time.time()

print "Finished in %.3f seconds" %(end - start)
"""


#bam_folder = sys.argv[1]
#run_locally = sys.argv[2]
#BAMS = open(bam_folders,'r')
#for bam_folder in BAMS:

parser = argparse.ArgumentParser('Test Parser',description='Work Around for now..')
parser.add_argument('bam_folder',help='The bam folder location')
parser.add_argument('-run_locally','--run_locally',help='Whether to run the luigi job locally',action='store_true',dest='run_locally')
args = parser.parse_args()

bam_folder = args.bam_folder
run_locally = args.run_locally


bam_folder = bam_folder.strip('\n')
sample = bam_folder.split('/')[-1]
bam_file = os.path.join(bam_folder,sample+'.realn.recal.bam')
gvcf = os.path.join(bam_folder,sample+'.g.vcf.gz')

#out_file = os.path.join(bam_folder,sample+'.genomecvg.bed')
cvg_folder = os.path.join(bam_folder,'cvg_binned')
gq_folder = os.path.join(bam_folder,'gq_binned')

print bam_file
if os.path.isfile(bam_file):
## Run Luigi

    luigi_cvg_cmd=(
        """/nfs/goldstein/software/python2.7.7/bin/luigi"""
        """ --module qc_module Binning --bam-file {0}"""
        """ --sample {1} --scratch {2} --mode coverage --block-size 10000 --poll-time 60 --worker-wait-interval 60 --no-tarball""".format(bam_file,sample,cvg_folder)
    )

    luigi_gq_cmd=(
        """/nfs/goldstein/software/python2.7.7/bin/luigi"""
        """ --module qc_module Binning --gvcf {0}"""
        """ --sample {1} --scratch {2} --mode gq --block-size 10000 --poll-time 60 --worker-wait-interval 60 --no-tarball""".format(gvcf,sample,gq_folder)
    )

    args = shlex.split(luigi_cvg_cmd)
    print args
    args = shlex.split(luigi_gq_cmd)
    print args
    os.system(luigi_gq_cmd)
    os.system(luigi_cvg_cmd)
#/nfs/goldstein/software/python2.7.7/bin/luigi --module qc_module CreateGenomeBed --bam /nfs/fastq16/ALIGNMENT/BUILD37/DRAGEN/EXOME/washei40496/washei40496.bam --sample washei40496 --scratch /nfs/fastq16/ALIGNMENT/BUILD37/DRAGEN/EXOME/washei40496
#BAMS.close()



