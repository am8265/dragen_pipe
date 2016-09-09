#!/usr/bin/env Python

import argparse
import gzip
import bz2
import os
from collections import defaultdict,OrderedDict

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
    """ Checks whether the alt allele is the string '<NON_REF>', if so return True else False
    
    Args : alt ; string ; the alternate allele column in the vcf file 
    Returns : bool  
    """

    return (alt == 'NON_REF')

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
    except ValueError as e: ## Need to look for stray cases where DP is not defined in the info field, it is (atleast in the cases i have encountered) 0 in such cases 
        return 0

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
        print Exception("Warning ! Caught a multiallelic flag through the filter-Check code/vcf")
        return [-1,-1]

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
    
def return_pedmap(gvcf_hashmap,coordinates):
    """ Create the genotypes and loci for the ped map files by processing the gvcf_hashmap at the input coordinates

    Args: gvcf_hashmap; dictionary of lists ; key : 'chrom_pos', val : [rsid,ref,alt,gt,dp,variant_type]
          coordinates; dictionary; key: 'chrom_pos', val: 0/1
    Returns: List, List of tuples ; the ped file genotypes, the map file info ; these can later be processed for outputing a ped file
    """
    
    temp_ped = []
    temp_map = []
        
    for key in coordinates:
        chrom,pos = key.split('_')                    
        if key in gvcf_hashmap: ## Site is present in gvcf file
            rsid,ref,alt,gt = gvcf_hashmap[key][0:4]
            if gvcf_hashmap[key][-1] == 'nonref': ## Non ref site
                if int(gvcf_hashmap[key][-2]) >= 10: ## Check dp
                    temp_ped.extend([ref,ref])
                    temp_map.append((chrom,rsid,'0',pos))
            else: ## A variant site
                if not is_multiallelic(gt) and is_valid_snv(ref,alt):
                    allele1,allele2 = get_genotypes(ref,alt,gt)
                    if allele1 == -1 and allele2 == -1: ## Stringent check,set to missing
                        temp_ped,temp_map = update_missing_variants(temp_ped,temp_map,(chrom,rsid,'0',pos))
                    temp_ped.extend(get_genotypes(ref,alt,gt))
                    temp_map.append((chrom,rsid,'0',pos))
                else: ## Consider the site to be missing
                    temp_ped,temp_map = update_missing_variants(temp_ped,temp_map,(chrom,rsid,'0',pos))
        else: ## Site is missing
            rsid = coordinates[key]
            temp_ped,temp_map = update_missing_variants(temp_ped,temp_map,(chrom,rsid,'0',pos))
            
    return [temp_ped,temp_map]
                                        
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

    contents = line.strip('\n').split('\t')
    chrom,pos,rsid,ref,alt = contents[0:5]
    info = contents[-2]
    info_contents = info.split(':')
    if len(info_contents) == 5: ## A non ref site
        end = contents[-3].split('=')[1]
        alt_type = 'nonref'
    else: ## A ref site
        end = pos
        alt_type = 'alt'
    val_field = contents[-1]
    gt = get_info_key(info,'GT',val_field)
    dp = get_info_key(info,'DP',val_field)
    interval = int(end) - int(pos) + 1

    if alt_type == 'alt':
        yield [chrom,pos,rsid,ref,alt,gt,dp,'alt']
    else:
        for _ in range(interval):
            yield [chrom,pos,rsid,ref,alt,gt,dp,'nonref']
            pos = str(int(pos) + 1)

def print_output(line,output_file):
    """ Print output line to file

    Args : line ; string ; the line(s) to print
    """
    
    with open(output_file,'w') as OUT:
        OUT.write(line+'\n')
           
def iterate_gvcf(gvcf):
    """ Iterate over the gvcf file and keep storing the genotype and marker loci 
    information
    
    Args : gvcf ; string ; path to the gvcf file 
    Returns : dict;
    """

    genotypes_ped = []
    locations_map = []

    gvcf_hashmap = {}
    
    with open_by_magic(gvcf) as GVCF:
        for line in GVCF:
            if line[0]!='#': ## Skip header
                gvcf_generator = parse_gvcf(line) ## If the site is multiallelic there could be more than one return 
                for element in gvcf_generator: 
                    chrom,pos,rsid,ref,alt,gt,dp,alt_type = element
                    ## Better to store only valid sites so as to reduce space complexity
                    if alt.find(',')!=-1: ## A multi-allelic site, use the first allele, this could be the alt,<NON_REF> case , but these are still valid since the GT field signfies these to have only 1 valid alt allele, so lets use them ! 
                        alt = alt.split(',')[0]               
                        gvcf_hashmap[chrom+'_'+pos] = [rsid,ref,alt,gt,dp,alt_type]
                
    return gvcf_hashmap

def parse_map_file(map_file):
    """ Function to parse a map file and return the ordered positions 
    Args : map_file ; string ; the map file
    Returns : list ; the coordinates in the map file
    """

    #coordinates_info = defaultdict(lambda: -1)
    coordinates = OrderedDict()
    with open(map_file,'r') as MAP:
        for line in MAP:
            contents = line.strip('\n').split('\t')
            pos = contents[-1]
            chrom = contents[0]
            rsid = contents[1]
            #coordinates[chrom+'_'+pos] = 0
            coordinates[chrom+'_'+pos] = rsid
    
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
    
    #genotypes_ped,locations_map = iterate_gvcf(gvcf,coordinates)
    gvcf_hashmap = iterate_gvcf(gvcf)
    print "Loaded gvcf hashmap"
    genotypes_ped,locations_map = return_pedmap(gvcf_hashmap,coordinates)
    
    
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
