#!/usr/bin/env Python

import argparse
import gzip
import bz2
import os
from collections import defaultdict,OrderedDict
import pysam
import warnings
import logging 
from pyfaidx import Fasta

########################################################################
###                                                                  ###
###              Script for creating ped,map file                    ### 
###              from a vcf and bam file                             ###
###              Author : Raghav                                     ###
########################################################################


def get_coverage_from_bam(BAM,chrom,pos):
    """ Get the coverage at a site, i.e. chromosome, position
    from a bam file

    BAM : pysam object ; the bam file opened by pysam
    chrom : str ; the chromosome 
    pos : int ; the position

    returns : int ; the coverage at the site
    """

    
    for pileupcolumn in BAM.pileup(chrom,pos-1,pos+1):
        if pileupcolumn.pos == pos-1: ## Match the exact site,
                ## since bam files are 0-based and vcfs are 1-based,the site for which
                ## we want the coverage will be at (pos-1) in the bam file 
            return pileupcolumn.n
            

def is_valid_snv(ref,alt):
    """ Checks to make sure that the snv is valid for a ped file
    
    Args : ref ; string ; the ref allele
           alt ; string ; the alt allele
    Returns : bool ; True if the conditions are met
    """

    if len(ref)!=1 or len(alt)!=1:
        return False
    if set(ref).issubset('ATCG') and set(alt).issubset('ATCG'):
        return True
    else:
        return False
    
def open_by_magic(filename):
    """ Taken from : http://stackoverflow.com/questions/18367511/how-do-i-automatically-handle-decompression-when-reading-a-file-in-python
    Uses the initial bytes of a file to detect the file compression.

    Args : filename ; string ; Path to the file
    Returns : file handle to the input file
    """
    magic_dict = {
            "\x1f\x8b\x08": gzip.open,
            "\x42\x5a\x68": bz2.BZ2File,
        }
    max_len = max(len(x) for x in magic_dict)
    
    with open(filename) as f:
        file_start = f.read(max_len)
    for magic, fn in magic_dict.items():
        if file_start.startswith(magic):
            return fn(filename)
    return open(filename,'r') ## Otherwise just a regular file
           
def is_multiallelic(gt):
    """ Returns a boolean based on whether a variant site is multiallelic
    This is deciphered using the GT field in the gvcf/vcf file

    Args : gt ; string ; e.g. 0/0,1/1,1/0,2/0,etc.etc.
    Returns : bool ; True if the site is multiallelic, False if not
    """
    
    if set(gt.split('/')).issubset(set(['0','1','2','3'])): ## There are cases when there is a '.' in the GT field
        return (max([int(allele) for allele in gt.split('/')]) > 1)
    else:
        return True

def get_info_key(info,key,field):
    """ The last 2 columns in a vcf/gvcf file are like : 
    GT:DP:GQ:MIN_DP:PL    0/0:2:6:20,6,65
    Here GT would be 0/0, DP would be 2 , GQ 6 and so on
    There is not set ordering required hence the values in the last column
    are ordered according to the names in the penultimate column. This function returns
    the position of the info key(say GT) in a given info field. So for the above e.g. if
    the key specified was GQ this function will return 2(the indexing is 0 based), this can 
    then be used to query the final column for the requisite key value

    Args : info ; string ; the penultimate column of a vcf/gvcf file 
           key ; string ; the specific info key for which we need the index
    Returns : int ; the index(0 based) of the key in the info column   
    """

    try:
        return field.split(':')[info.split(':').index(key)]
    except ValueError as e: ## Need to look for stray cases where DP is not defined in the info field, it is 0 in such cases  (atleast in what I have encountered) 
        return 0

def get_genotypes_bam(ref,dp):
    """ If the coverage is < 10 the site is treated as missing
    otherwise we treat it as a homozygous reference
    
    ref : str ; the reference allele 
    dp : str ; the depth at the site    

    return : list ; [allele1,allele2]
    """

    if dp >= 10:
        return [ref,ref]
    else:
        return ['0','0']
    
def get_genotypes(ref,alt,gt):
    """
    Uses the GT field of a vcf files and returns the appropriate genotypes for the 2 alleles
    Note : I am assuming the ref allele to be the allele1 in the plink file 

    Args: ref; string ; the reference allele genotype
          alt; string ; the alternate allele genotype
          gt; the GT field in the vcf file specifying whether the variant is hom_ref,hom_alt,etc.
    """
    if gt == '0/0':
        return [ref,ref]
    elif gt == '1/1':
        return [alt,alt]
    elif gt == '1/0' or gt == '0/1':
        return [ref,alt]
    else: ## Multi allelic site
        raise Exception("Caught a multiallelic flag through the filter-Check code/vcf")

def update_missing_variants(ped_list,map_list,map_info):
    """ Helper function for updating missing entries
    Args: ped_list ; list ; contains output ped info
          map_list ; list of tuples ; contains output map info
          map_info ; tuple ; contains info for the map file
    Returns: list,list of tuples; update and return input data structures 
    """
    ped_list.extend(['0','0'])
    map_list.append(map_info)
    return [ped_list,map_list]
    
def return_pedmap(vcf_hashmap,coordinates,bam,genome):
    """ Create the genotypes and loci for the ped map files by processing 
    the vcf_hashmap at the input coordinates

    Args: vcf_hashmap: dictionary of lists ; key : 'chrom_pos', val : [rsid,ref,alt,gt,dp]
          coordinates: dictionary; key: 'chrom_pos', val: 0/1
          bam : pysam object ; the bam file opened by pysam
          genome : pyfasta object ; the reference fasta 

    Returns: List, List of tuples ; the ped file genotypes, the map file info ; these can later be processed for outputing a ped file
    """
    
    temp_ped = []
    temp_map = []
    temp_geno = [] ## For .geno files (for pca and ethnicity)
    i = 0
    for key in coordinates:
        chrom,pos = key.split('_')
        id_var = coordinates[chrom+'_'+pos]
        if chrom+'_'+pos in vcf_hashmap: ## Site is present in vcf file
            ref,alt,gt = vcf_hashmap[key][0:3]
            if not is_multiallelic(gt) and is_valid_snv(ref,alt):
                temp_ped.extend(get_genotypes(ref,alt,gt))
                temp_map.append((chrom,id_var,'0',pos))
                if gt == '1/1':
                    temp_geno.append('2')
                else:
                    temp_geno.append('1')
            else:
                logging.warning("Warning : Multi allelic/Indel/MNP site : ref : {0}; alt : {1} ; gt : {2}".format(ref,alt,gt))
                temp_ped,temp_map = update_missing_variants(temp_ped,temp_map,(chrom,id_var,'0',pos))
                i+=1
        else:
            bam_coverage = get_coverage_from_bam(bam,chrom,int(pos))
            if bam_coverage < 10:
                i+=1
                temp_geno.append('9')
            else:
                temp_geno.append('0')
            ref = get_ref_allele(genome,chrom,int(pos))
            temp_ped.extend(get_genotypes_bam(ref,bam_coverage))
            temp_map.append((chrom,id_var,'0',pos))

    logging.info("{0} site(s) were missing".format(i))
    return [temp_ped,temp_map,temp_geno]


def get_ref_allele(genome,chrom,pos):
    """ Get the reference site from a pyfasta object
    
    genome : a pyfasta object ; the refernce genome
    chrom : str ; the chromosome 
    pos : int ; the position
    
    return : str ; the reference allele 
    """

    return str(genome[chrom][pos-1:pos]) ## Again the fasta file will be 0-based,
                                         ## so pos will correspond to pos-1 in the fasta

def parse_vcf(line):
    """ Parse a vcf line to get : 
    Chromosome , Position, Ref, Alt, GT
    Downstream function will ignore all sites which are 
    indels. For multi allelic sites, only the first allel
    is kept 

    Args ; line ; string ; The gvcf line 
    Returns ; list ; [chromosome, position, reference allele, alternate allele, genotype(GT)]
    """

    contents = line.strip('\n').split('\t')
    chrom,pos,rsid,ref,alt = contents[0:5]
    info = contents[-2]
    info_contents = info.split(':')    
    val_field = contents[-1]
    gt = get_info_key(info,'GT',val_field)
    dp = get_info_key(info,'DP',val_field)
    return [chrom,pos,ref,alt,gt,dp]

def print_output(line,output_file):
    """ Print output line to file

    Args : line ; string ; the line(s) to print
    """
    
    with open(output_file,'w') as OUT:
        OUT.write(line+'\n')
           
def iterate_vcf(vcf):
    """ Iterate over the vcf file and keep storing the genotype and marker loci 
    information
    
    Args : vcf ; string ; path to the gvcf file 
    Returns : dict;
    """

    vcf_hashmap = {}
    
    with open_by_magic(vcf) as VCF:
        for line in VCF:
            if line[0]!='#': ## Skip header
                chrom,pos,ref,alt,gt,dp = parse_vcf(line)
                if alt.find(',')!=-1: ## A multi allelic site,might need to look at GT more carefully for this
                    alt = alt.split(',')[0]
                    gt = gt.split(',')[0]
                vcf_hashmap[chrom+'_'+pos] = [ref,alt,gt,dp]

    return vcf_hashmap

def parse_map_file(map_file):
    """ Function to parse a map file and return the ordered positions 
    Args : map_file ; string ; the map file
    Returns : list ; the coordinates in the map file
    """

    coordinates = OrderedDict()
    with open(map_file,'r') as MAP:
        for line in MAP:
            contents = line.strip('\n').split('\t')
            pos = contents[3]
            chrom = contents[0]
            id_var = contents[1]
            coordinates[chrom+'_'+pos] = id_var
    
    return coordinates
            
def main(args):
    """ Main Function
    """

    ## Unpack the arguments
    sample_id = args.sample_id
    sample_id = sample_id + '_' + args.pseudo_prepid
    if args.family_id == 'dummy': ## User chose not to define a family id use the sample id as the family id
        family_id = sample_id
    else:
        family_id = args.family_id
    vcf = args.vcf
    markers = args.markers
    output_stem = args.stem
    sex = args.sex
    pheno = args.pheno
    if sex == 'M':
        sex = '1'
    else:
        sex = '2'
    bam = pysam.AlignmentFile(args.bam,"rb")    
    ## Create output directory if it doesnt exist
    out_dir = args.out
    if out_dir != './':
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)           
    ## Set logging
    logging.basicConfig(filename=args.log_file,level=logging.DEBUG,
                        filemode='w')
    ## The input markers
    coordinates = parse_map_file(markers)    
    vcf_hashmap = iterate_vcf(vcf)
    logging.info("Loaded vcf hashmap")
    ## Load the reference genome
    genome = Fasta(args.ref_fasta)
    logging.info("Loaded reference genome")
    genotypes_ped,locations_map,genotypes_geno = return_pedmap(vcf_hashmap,coordinates,bam,genome)
    
    temp = [family_id,sample_id,'0','0',sex,pheno]
    temp.extend(genotypes_ped)
    output_ped_line = ' '.join(temp)
    ## Unpack the locational information for the map file
    temp = []
    for loc in locations_map:
        line = '\t'.join(str(e) for e in loc)
        temp.append(line)
    output_map_line = '\n'.join(temp)

    logging.info("Writing the ped file")
    ## Print output files
    output_ped = os.path.join(out_dir,output_stem+'.ped')
    output_map = os.path.join(out_dir,output_stem+'.map')
    print_output(output_ped_line,output_ped)
    print_output(output_map_line,output_map)

    logging.info("Writing the geno file")
    output_geno_line = '\n'.join(genotypes_geno)
    output_geno = os.path.join(out_dir,output_stem+'.geno')
    print_output(output_geno_line,output_geno)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser('Create a ped file from a vcf file and a bam file',description='Takes as input : \n1. A valid gvcf file\n2. Marker loci\nOutput : A ped file containing genotype information at the appropriate marker loci')
    parser.add_argument('--sample_id','-sample_id',dest='sample_id',help='CHGVID',required=True)
    parser.add_argument('--family_id','-family_id',dest='family_id',required=False,help='The family id (if not given then we use the sample_id as the family id)',default='dummy')
    parser.add_argument('--vcf','-vcf',dest='vcf',help='The input vcf',required=True)
    parser.add_argument('--markers','-markers',dest='markers',help='The input genetic loci, the file should adhere to the plink map format : http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map\nIn short there should be exactly 4 columns :\n1.Chromosome\n2.Snp identifier\n3.Genetic distance\n4.Base-pair position (bp units)\nFor this program it is enough if you have the 1st and 4th column filed in correctly, just fill in the other 2 columns with dummy data if you do not have that information',required=True)
    parser.add_argument('--sex','-sex',dest='sex',help='M/F',required = False,default='M',choices=['M','F'])
    parser.add_argument('--pheno','-pheno',dest='pheno',help='1/2',required = False,default='1',choices=['1','2'])
    parser.add_argument('--bam','-bam',dest='bam',help="Path to the bam file",required = True)
    ## Note I am making prepid and seqtype as mandatory ,
    ## since at IGM these attributes will ensure uniqueness of sample
    ## which is essential when you do downstream analysis with ped files
    ## The sample name will be a concatenation of sample_id, seqtype and prepid, for e.g. chgvid_seqtype_prepid
    parser.add_argument('--seqtype','-seq',dest='seq',help='The type of sequencing, i.e. Exome or Genome', required = True)
    parser.add_argument('--pseudo_prepid','-pseudo_prepid',dest='pseudo_prepid',help='The prepid for your sample',required=True)
    parser.add_argument('--outdir','-out',dest='out',help='Output dir, can give full path, output directories will be created if it doesnt exist',required=True)
    parser.add_argument('--stem','-stem',dest='stem',help='The stem for the output file, a .ped and .map will be appended to this for the output files',required=True)
    parser.add_argument('--ref','-ref',dest='ref_fasta',help='The reference genome',required=False,default='/scratch/HS_Build37/BWA_INDEX_hs37d5/hs37d5.fa')
    parser.add_argument('--logfile','-log',dest='log_file',help='The full path to the logfile',required=False,default='log.txt')
    args = parser.parse_args()
    main(args)
