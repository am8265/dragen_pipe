#!/usr/bin/env Python

import argparse
import gzip
import bz2
import os
from collections import defaultdict

########################################################################
###                                                                  ###
###              Script for creating ped,map file                    ### 
###              from a gvcf file                                    ###
###              Author : Raghav                                     ###
########################################################################


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

def is_nonref_alt(alt):
    """ Checks whether the alt allele is the string <NON_REF>, if so return True else False
    
    Args : alt ; string ; the alternate allele column in the vcf file 
    Returns : bool  
    """

    return (alt == 'NON_REF')
    
def parse_gvcf(line):
    """ Parse a gvcf line to get : 
    Chromosome , Position, Ref, Alt, GT
    Will ignore all sites which are multiallelic and have only a <NON_REF> as the alternate
    allele, there are sites where there is a valid alt allele and a <NON_REF>
    allele as well, these will not be treated as multi allelic sites since the 
    GT field will be either 0/0, 1/1 or 1/0, as a general rule if a site is valid
    i.e. not multiallelic and does not have only a <NON_REF> alt allele, then we will use the first
    allele in the alt column, for now further checks to DP/GQ are not being enforced 

    Args ; line ; string ; The gvcf line 
    Returns ; list ; [chromosome, position, reference allele, alternate allele, genotype(GT)]
    """

    chrom,rsid,pos,ref,alt = line.strip('\n').split('\t')[0:5]
    gt = line.strip('\n').split('\t')[-1].split(':')[0] ## Maybe generalize this later based on the string GT in the penultimate column
    return [chrom,rsid,pos,ref,alt,gt]

def print_output(line,output_file):
    """ Print output line to file

    Args : line ; string ; the line(s) to print
    """

    with open(output_file,'w') as OUT:
        OUT.write(line+'\n')
           
def iterate_gvcf(gvcf,positions):
    """ Iterate over the gvcf file and keep storing the genotype and marker loci 
    information
    
    Args : gvcf ; string ; path to the gvcf file 
    
    Returns : (list,tuple) ; the genotypes and the marker loci
    """

    genotypes_ped = []
    locations_map = []
    
    with open_by_magic(gvcf) as GVCF:
        for line in GVCF:
            if line[0]!='#': ## Skip header
                chrom,pos,rsid,ref,alt,gt = parse_gvcf(line)
                if alt.find(',')!=-1: ## Potentially multi-allelic, use the first allele
                    alt = alt.split(',')[0]
                ## Check if the site is valid
                if (not is_nonref_alt(alt) and  not is_multiallelic(gt) and is_valid_snv(ref,alt) and positions[chrom+'_'+pos] == 0):
                    if gt == '0/0':
                        genotypes_ped.extend((ref,ref))
                        locations_map.append((chrom,rsid,'0',pos))
                    elif gt == '1/0':
                        genotypes_ped.extend((ref,alt))
                        locations_map.append((chrom,rsid,'0',pos))
                    elif gt == '1/1':
                        genotypes_ped.extend((alt,alt))
                        locations_map.append((chrom,rsid,'0',pos))

    return (genotypes_ped,locations_map)

def parse_map_file(map_file):
    """ Function to parse a map file and return the positions 
    Args : map_file ; string ; the map file
    Returns : list ; the coordinates in the map file
    """

    coordinates = defaultdict(lambda: -1)
    with open(map_file,'r') as MAP:
        for line in MAP:
            pos = line.strip('\n').split('\t')[-1]
            chrom = line.strip('\n').split('\t')[0]
            coordinates[chrom+'_'+pos] = 0 
    
    return coordinates
            
def main(args):
    """ Main Function
    """

    ## Unpack the arguments 
    sample_id = args.sample_id
    if args.family_id == 'dummy': ## User chose not to define a family id use the sample id as the family id
        family_id = sample_id
    gvcf = args.gvcf
    markers = args.markers
    output_stem = args.output
    sex = args.sex
    pheno = args.pheno
    if sex == 'M':
        sex = '1'
    else:
        sex = '2'

    ## Create output directory if it doesnt exist
    out_dir = os.path.dirname(output_stem)
    if out_dir != '':
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)

    ## The input markers
    coordinates = parse_map_file(markers)
    
    genotypes_ped,locations_map = iterate_gvcf(gvcf,coordinates)
   
    temp = [sample_id,family_id,'0','0',sex,pheno]
    temp.extend(genotypes_ped)
    output_ped_line = ' '.join(temp)
    ## Unpack the locational information for the map file
    temp = []
    for loc in locations_map:
        line = '\t'.join(str(e) for e in loc)
        temp.append(line)
    output_map_line = '\n'.join(temp)

    ## Print output files
    print_output(output_ped_line,output_stem+'.ped')
    print_output(output_map_line,output_stem+'.map')
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser('Create a ped file from a gvcf file',description='Takes as input : \n1. A valid gvcf file\n2. Marker loci\nOutput : A ped file containing genotype information at the appropriate marker loci')
    parser.add_argument('-sample_id','--sample_id',dest='sample_id',help='CHGVID',required=True)
    parser.add_argument('-family_id','--family_id',dest='family_id',required=False,help='The family id (if not given then we use the sample_id as the family id)',default='dummy')
    parser.add_argument('-gvcf','--gvcf',dest='gvcf',help='The input gvcf',required=True)
    parser.add_argument('--markers','-markers',dest='markers',help='The input genetic loci, the file should adhere to the plink map format : http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map\nIn short there should be exactly 4 columns :\n1.Chromosome\n2.Snp identifier\n3.Genetic distance\n4.Base-pair position (bp units)\nFor this program it is enough if you have the 1st and 4th column filed in correctly, just fill in the other 2 columns with dummy data if you do not have that information',required=True)
    parser.add_argument('--sex','-sex',dest='sex',help='M/F',required = False,default='M',choices=['M','F'])
    parser.add_argument('--pheno','-pheno',dest='pheno',help='1/2',required = False,default='1',choices=['1','2'])
    parser.add_argument('-output','--output',dest='output',help='Prefix for output ped and map files, can give full path, output directories will be created if it doesnt exist',required=True)
    args = parser.parse_args()
    main(args)
