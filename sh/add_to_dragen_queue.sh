[config]
#Setup
gatk=/nfs/goldstein/software/GATK-3.5.0
java=/nfs/goldstein/goldsteinlab/software/java/jdk1.7.0_03/bin/java
max_mem=15
ref=/nfs/goldsteindata/refDB/HS_Build37/BWA_INDEX_hs37d5_BWAmem/hs37d5.fa
#scratch=/nfs/seqscratch09/ALIGNMENT/DRAGEN/EXOME/%(sample_name)s
scratch=/nfs/seqscratch09/ALIGNMENT/DRAGEN
base_directory=/nfs/fastq_temp/ALIGNMENT/BUILD37/DRAGEN/EXOME
#sample_name=fetal0007M
interval_list=%(scratch)s/%(sample_name)s.interval_list
recal_table=%(scratch)s/%(sample_name)s.recal_data.table

#BAMs
bam=%(base_directory)s/%(sample_name)s/%(sample_name)s.bam
realn_bam=%(base_directory)s/%(sample_name)s/%(sample_name)s.realn.bam
#realn_recal_bam

#VQSR
vcf=%(scratch)s/%(sample_name)s.raw_vcf
snp_vcf=%(scratch)s/%(sample_name)s.snp_vcf
#snp_recal
#snp_tranches
#snp_rscript
#snp_filtered

indel_vcf=%(scratch)s/%(sample_name)s.indel_vcf
#indel_rscript
#indel_recal
#indel_tranches
#indel_rscript
#indel_filtered
#indel_filter_matched

#resources
hapman=/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/hapmap_3.3.b37.vcf 
omni=/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/1000G_omni2.5.b37.vcf 
g1000=/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/1000G_phase1.indels.b37.vcf 
Mills1000g=/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/Mills_and_1000G_gold_standard.indels.b37.vcf
dbSNP=/nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/dbsnp_138.b37.vcf

## Post GATK Variables
#java=/nfs/goldstein/software/jdk1.8.0_05/bin/java
picard=/nfs/goldstein/software/picard-tools-2.2.1/picard.jar
#reference_file=/nfs/goldsteindata/refDB/HS_Build37/BWA_INDEX_hs37d5/hs37d5.fa
seqdict_file=/nfs/goldsteindata/refDB/HS_Build37/BWA_INDEX_hs37d5/hs37d5.dict
bed_file=/home/rp2801/test_folder/coverage_statistics/ccds_r14.bed
target_file=/home/rp2801/test_folder/coverage_statistics/CCDS.genes.unique.regions.list
create_targetfile=False
bait_file=/home/rp2801/test_folder/coverage_statistics/SeqCap_EZ_Exome_v3_capture.list
wgsinterval=False
#scratch=/nfs/seqscratch11/rp2801/
#output_dir=./
#bam_file=/nfs/seqscratch11/filtered.Schiz6113363.production.bam
#sample_name=Schiz6113363
seq_type=Genome
python_path=/nfs/goldstein/software/python2.7.7/bin/python2.7
relatedness_refs=/nfs/goldstein/goldsteinlab/Bioinformatics/scripts/refs
sampleped_loc=/home/rp2801/git/bioinfo_pipe/Python/sampleped_fixed.py
relatedness_markers=/nfs/seqscratch09/RELATEDNESS_QC/test.map
bedtools_loc=/nfs/goldstein/software/bedtools-2.25.0/bin/bedtools
pypy_loc=/nfs/goldstein/software/pypy-4.0.1-linux_x86_64-portable/bin/pypy
coverage_binner_loc=/home/rp2801/git/dragen_pipe/dragen_coverage_binner.py