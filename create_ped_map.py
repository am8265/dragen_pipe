#!/usr/bin/env Python

import argparse
import gzip
import bz2

########################################################################
###                                                                  ###
###              Script for creating a ped file                      ### 
###              from a gvcf file                                    ###
###              Author : Raghav                                     ###
########################################################################


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

    return (max([int(allele) for allele in gt.split('/')]) <= 1)

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
        OUT.write(line)
           
def iterate_gvcf(gvcf):
    """ Iterate over the gvcf file and keep storing the genotype and marker loci 
    information
    
    Args : gvcf ; string ; path to the gvcf file 
    
    Returns : (list,tuple) ; the genotypes and the marker loci
    """

    genotypes_ped = []
    locations_map = []
    
    with open_by_magic(gvcf) as GVCF:
        for line in GVCF:
            chrom,rsid,pos,ref,alt,gt = parse_gvcf(line)
            ## Check if the site is valid
            if (not is_nonref_alt(alt) and  not is_multiallelic(gt)):
                if gt == '0/0':
                    genotypes_ped.append(ref,ref)
                    locations_map.append((chrom,rsid,pos))
                elif gt == '1/0':
                    genotypes_ped.append(ref,alt)
                    locations_map.append((chrom,rsid,pos))
                elif gt == '1/1':
                    genotypes_ped.append(alt,alt)
                    locations_map.append((chrom,rsid,pos))

    return (genotypes_ped,locations_map)

def main():
    """ Main Function
    """

    ## Unpack the arguments 
    sample_id = args.vcf
    if args.family_id == 'dummy': ## User chose not to define a family id use the sample id as the family id
        family_id = sample_id
    gvcf = args.gvcf
    markers = args.markers
    output = args.output

    ## Define keyword arguments for parsing to functions 
    kwargs = {'sample_id':sample_id,'gvcf':gvcf,'markers':markers,'output':output}

    genotypes_ped,locations_map = iterate_gvcf(gvcf)

    ped_line = ['sample_id'
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser('Create a ped file from a gvcf file',description='Takes as input : \n1. A valid gvcf file\n2. Marker loci\nOutput : A ped file containing genotype information at the appropriate marker loci')
    parser.add_argument('-sample_id','--sample_id',dest='sample_id',help='CHGVID',required=True)
    parser.add_argument('-family_id','--family_id',dest='family_id',required=False,help='The family id (if not given then we use the sample_id as the family id)',default='dummy')
    parser.add_argument('-gvcf','--gvcf',dest='gvcf',help='The input gvcf',required=True)
    parser.add_argument('--markers','-markers',dest='markers',help='The input genetic loci, the file should adhere to the plink map format : http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map\nIn short there should be exactly 4 columns :\n1.Chromosome\n2.Snp identifier\n3.Genetic distance\n4.Base-pair position (bp units)\nFor this program it is enough if you have the 1st and 4th column filed in correctly, just fill in the other 2 columns with dummy data if you do not have that information',required=True)
    parser.add_argument('--sex','-sex',dest='sex',help='M/F',required = False,default='M')
    parser.add_argument('-output','--output',dest='output',help='The destination for the output ped file',required=True)
    args = parser.parse_args()
    main(args)
