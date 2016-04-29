#!/nfs/goldstein/software/python2.7.7/bin/python
# -*- coding: utf-8 -*-
"""
Create dragen configuration file for gvcf creation from bam alignment
"""

from datetime import datetime
def create_config(sample):

    raw_config="""
#================================================================================
# Dragen 1.7.1 Configuration File
#================================================================================
#
# Author:	Joshua Bridgers (jb3816)
# Date:		02/04/2016
# Desc:         Configuration file for Dragen 1.7.1 for gVCF production from
#               a bam file
#
#================================================================================

#================================================================================
# SAMPLE SETUP
#================================================================================

ref-dir = /staging/REF/b37_decoy/
intermediate-results-dir = /staging/tmp
bam-input = {bamFile}

#================================================================================
# OUTPUT
#================================================================================

output-file-prefix = {IGMID}
output-directory = /nfs/fastq16/ALIGNMENT/BUILD37/DRAGEN/{uppercaseSeqtype}/{IGMID}/

#================================================================================
# VARIANT CALLING
#================================================================================

enable-variant-caller = true
vc-sample-name = {IGMID}
vc-reference = /staging/REF/b37_decoy/hs37d5.fa  	# same as --ref-dir parameter
vc-min-call-qual = 20 						# Default 30
vc-min-read-qual = 20 						# Min MAPQ. #Default 20
vc-min-base-qual = 10 						# Default 10
vc-min-reads-per-start-pos = 5 					# default 5; don’t downsample below 5 at a given starting position
vc-max-reads-per-region = 1000 					# Maximum number of reads per region for downsampling.  Default 1000
vc-target-coverage = 2000 					# The target coverage for downsampling.  Default 2000
dbsnp = /staging/REF/dbSNP_b144_GRCh37p13/All_20150605.vcf
vc-emit-ref-confidence = GVCF				#Activates GVCF mode

#================================================================================
# NOTES
#================================================================================
# Defaults taken from:
# Dragen User Defaults Configuration File - Revision 01.002.060.01.01.14.12343
# Adapter masking is performed during BCL
# No hard filtering applied because we handle that from the ATAV side
#
# http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/ #Patches to GRCh37 can be found here
# http://www.ncbi.nlm.nih.gov/variation/docs/multi_assembly_support/ #Multi Assembly Support
# http://exac.broadinstitute.org/faq
#
# vcf hard filter settings (Same as GATK)
# vc-hard-filter = DRAGENHardSNP:snp: QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0;DRAGENHardINDEL:indel: QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0

#================================================================================
# FIN
#================================================================================
"""

    outfile = output_config(sample,raw_config)
    return outfile

def output_config(sample,raw_config):
    final_config = raw_config.format(**sample.metadata)
    out_file = "{conf_file}".format(**sample.metadata)

    print out_file
    out = open(out_file,'w')
    out.write(final_config)
    out.close()
    return out_file

