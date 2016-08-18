import sys
import subprocess
from utilities import is_gzipped,fh
import argparse
import os


########################################################################
###                                                                  ###
###              Script for binning :                                ###
###      1. GQ values from a gvcf file                               ###
###      2. Coverage values from bedtools genomecvg file             ###
###              Author : Raghav                                     ###
########################################################################


def gq_to_letter(gq):
    
    """ Returns Letter encoding for gq value
        gq => numerical gq value
    """

    gq = int(gq)

    if(gq < 5):
        return 'a'
    elif(gq < 15):
        return 'b'
    elif(gq < 20):
        return 'c'
    elif(gq < 60):
        return 'd'
    else:
        return 'e'

def coverage_to_letter(coverage):
    
    """ Returns Letter encoding for coverage value
        coverage => numerical coverage value
    """

    coverage = int(coverage)
    if(coverage < 3):
        return 'a'
    elif(coverage <= 9):
        return 'b'
    elif(coverage <= 19):
        return 'c'
    elif(coverage <= 200):
        return 'd'
    else:
        return 'e'

    
def bin_mode_to_letter(mode_val, mode_len):
    
    """ Return the printable form for the alpha-numeric gq/coverage bins
        i.e. mode_val = ['a','b','c'] mode_len = [200,300,500]
        200a,300b,500c is returned
        mode_val => a list of letter encoded gq values
        mode_len => a list of integers corresponding to the lengths of the
                        above letter encoded gq values
    """

    return_mode_string_list = []
    for c,l in zip(mode_val, mode_len):
        l = radix36_encoding(l)
        return_mode_string_list.append(l+c) 

    return ''.join(return_mode_string_list)


def update_mode_bins(n, mode, mode_len, mode_val):
    
    """ Update the coded gq/coverage values according to the letter encodings
        i.e. if n = 100 , mode = 'a' , mode_val = ['b','c','a'] and 
        mode_len = [100,200,200], the program will update the gq_len
        to [100,200,300] and mode_val to ['b','c','a']
        n => numerical value of gq/coverage
        mode => letter code for gq/coverage bin
        mode_len => list of numerical gq values seen in the current block
        mode_val => list of letter codes for gq/coverage bins seen in the current block
    """

    if(not mode_len):  # If list is empty
        mode_len.append(n)
        mode_val.append(mode)

    else:
        if(mode == mode_val[-1]):
            mode_len[-1] = mode_len[-1] + n
            
        else:
            mode_len.append(n)
            mode_val.append(mode)

    return [mode_len, mode_val]

def radix36_encoding(number):
    
    """ 
    number => numerical gq value in decimal
    returns base36 encoded value for the above
    Refer to : http://stackoverflow.com/questions/1181919/python-base-36-encoding
    """ 

    alphabet = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    if not isinstance(number, (int, long)):
        raise TypeError('number must be an integer')

    base36 = ''
    sign = ''

    if number < 0:
        sign = '-'
        number = -number

    if 0 <= number < len(alphabet):
        return sign + alphabet[number]

    while number != 0:
        number,i = divmod(number,len(alphabet))
        base36 = alphabet[i] + base36

    return str(sign + base36)
      

def radix36decode(number):
    """Decode the radix36 number to decimal
    """
    return int(number,36)

def is_gq_present(line):
    """
    Returns a boolean depending on whether GQ 
    is present in a gvcf line
    """

    if line.find('GQ') == -1: ## No GQ info is available
        return False

    else:
        return True


def parse_gvcf_gq(line):
    """Parses a gvcf line to return chromosome,start,stop and gq     
    """

    
    found = False
  
    contents = line.split('\t')
    check_site = contents[-2]
    gq_field = contents[-1]
    chromosome = contents[0]
    start = int(contents[1])

    temp_check = check_site.split(':')

    if len(temp_check) == 5: ## Non Variant site
        gq = gq_field.split(':')[2]
        end = int(contents[7].strip('END='))

    else: ## Variant site
        for i in range(0,len(temp_check)):
            if temp_check[i] == 'GQ':
                found = True
                gq_id = i
                break

        if found == True:
            gq = gq_field.split(':')[gq_id]
        else:
            gq = 0
            #print line  ## For debugging 
        end = start

    

    return [chromosome,start,end,gq]


def is_empty_modelist(mode_len):
    """
    Check whether the mode_bin length list is empty
    returns a boolean
    """
    if mode_len:
        return False 
    else:
        return True

def output_line(mode_str,exclusion_str,sample_id,out_start,block,bin_out):
    """
    Output gq_str to a file as long as it is not the
    exclusion string
    """

    if mode_str != exclusion_str:
        bin_out.write("%s\t%s\t%s\n" % (
                sample_id,str(out_start/block),mode_str))

        
def chromosome_check(prev_chromosome,chromosome):
    """
    Check whether the chromsome has changed 
    returns True/False
    """

    return prev_chromosome == chromosome 
                 
    
def reinitialize_variables(out_start,block):
    """
    Reinitializes key variables after a block has been spanned
    """

    mode_len = []
    mode_val = []
    out_start = out_start + block
    current_block = block
    
    return (out_start,current_block,mode_len,mode_val)


def parse_genomecvg(line):
    """
    line : tsv bed file line
    return : [chrom,start,end,coverage]
    return types : [str,int,int,int]
    """

    
    c = line.split('\t')
    return [c[0],int(c[1]),int(c[2]),int(c[3])]


def main(args):
    """ Main Function """

    
    block = args.B
    input_file = args.inputfile
    sample_id = args.sample_id
    out_dir = args.output_dir

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)


    IN = fh(input_file)
    
    human_chromosomes = range(1,23)
    human_chromosomes = [str(x) for x in human_chromosomes]
    human_chromosomes.extend(['X','Y','MT'])
    prev_chromosome = '1'
    exclusion_str = radix36_encoding(block) + 'a'
    flag = 0


    mode_len = []
    mode_val = []
    
    
    for line in IN:
        line = line.strip('\n')
        if args.mode == 'gq':
            if line[0] == '#':
                continue
            chromosome,start,end,mode = parse_gvcf_gq(line)
        if args.mode == 'coverage':
            chromosome,start,end,mode = parse_genomecvg(line)
        
        
        if chromosome not in human_chromosomes:
            continue ## Skip scaffolds
        
        if chromosome_check(prev_chromosome,chromosome) == False:
            if is_empty_modelist(mode_len) == False: ## Account for remaining bins if present
                if block - sum(mode_len) != 0:  ## Make sure the block is filled
                    update_mode_bins(
                        block - sum(mode_len),'a',mode_len,mode_val)
                
                    if args.mode == 'gq':
                        mode_str = bin_mode_to_letter(
                            mode_val,mode_len)
                    else: ## Not checking further since the choices are restriced in argparse
                        mode_str = bin_mode_to_letter(
                            mode_val,mode_len)
                                           
                    output_line(mode_str,exclusion_str,sample_id,out_start,block,bin_out)
                    out_start,current_block,gq_len,gq_val = reinitialize_variables(out_start,block)

            bin_out.close()
            
            bin_out = open(os.path.join(out_dir,sample_id+'_%s_binned_'%args.mode+str(block)+'_chr%s.txt'%chromosome),'w')
            
                
            flag = 0 
                

        if flag == 0:
            bin_out = open(os.path.join(out_dir,sample_id+'_%s_binned_'%args.mode+str(block)+'_chr%s.txt'%chromosome),'w')
            out_start = start - start % block
            current_block = block
            mode_len = []
            mode_val = []
            
        current_interval = int(end) - int(start) + 1     
    


        if args.mode == 'gq':
            mode_letter = gq_to_letter(mode)
        else:
            mode_letter = coverage_to_letter(mode)


        if current_block - current_interval > 0: ## Block spans the interval

            # Update block size
            current_block = current_block - current_interval
            n = current_interval
            
            mode_len,mode_val = update_mode_bins(
                n,mode_letter,mode_len,mode_val)
        
                

        elif current_block - current_interval == 0:
            n = current_interval
            mode_len,mode_val = update_mode_bins(
                n,mode_letter,mode_len,mode_val)
            if args.mode == 'gq':
                mode_str = bin_mode_to_letter(
                    mode_val,mode_len)
            else:
                mode_str = bin_mode_to_letter(
                    mode_val,mode_len)

            ## print output to file
            output_line(mode_str,exclusion_str,sample_id,out_start,block,bin_out)
        
            ## Reinitialization
            out_start,current_block,mode_len,mode_val = reinitialize_variables(
                out_start,block)
           
                
        elif current_block - current_interval < 0:
            n = current_block
            mode_len,mode_val = update_mode_bins(
                n,mode_letter,mode_len,mode_val)

            
            mode_str = bin_mode_to_letter(
                mode_val,mode_len)
            
            mode_str = bin_mode_to_letter(
                mode_val,mode_len)

            output_line(mode_str,exclusion_str,sample_id,out_start,block,
                        bin_out)

            ## Reinitialization
            out_start,current_block,mode_len,mode_val = reinitialize_variables(
                out_start,block)            
            remaining_interval = current_interval - current_block
            num_blocks = remaining_interval/block
            
            
            if(current_block - remaining_interval <= 0):
                remaining_interval = remaining_interval % current_block
                while(num_blocks > 0):
                    n = current_block
                    mode_len, mode_val = update_mode_bins(
                        n, mode_letter, mode_len, mode_val)

                    if args.mode == 'gq':
                        mode_str = bin_mode_to_letter(
                            mode_val, mode_len)
                    else:
                        mode_str = bin_mode_to_letter(
                            mode_val, mode_len)

                    output_line(mode_str,exclusion_str,sample_id,out_start,block,
                        bin_out)
                    
                    ## Reinitialization
                    out_start,current_block,mode_len,mode_val = reinitialize_variables(
                        out_start,block)
                    num_blocks = num_blocks - 1
                    

                if(remaining_interval > 0):
                    n = remaining_interval
                    mode_len, mode_val = update_mode_bins(
                        n, mode_letter, mode_len, mode_val)
                    # Update block size
                    current_block = current_block - remaining_interval

                else:
                    n = remaining_interval
                    mode_len, mode_val = update_mode_bins(
                        n, mode_letter, mode_len, mode_val)
                    # Update block size
                    current_block = current_block - remaining_interval
                  
        flag = 1
        # Store previous values
        prev_start = start
        prev_stop = end
        prev_chromosome = chromosome
        
        
    if mode_len : ## Account for remaining bins if present
        # Add a's if the block is not fully spanned
        sum_length = sum(mode_len)
        if block - sum_length != 0:
            mode_len,mode_val = update_mode_bins(block - sum_length, 'a',
                                           mode_len,mode_val)


        
        mode_str = bin_mode_to_letter(mode_val,mode_len)
        output_line(mode_str,exclusion_str,sample_id,out_start,block,
                        bin_out)
            
            
    IN.close()
    bin_out.close()

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Coverage/GQ Binning Script',description="Takes a genome coverage bed file produced by bedtools genomecov function or a GVCF file from GATK and creates upto 24 output files. Each file contains the base 36 encoded binned coverage/gq strings over a specified block interval for a particular chromosome. The script can be run with pypy, a JIT compiler to improve the run time by upto 10X (tested on pypy-2.2.1, though can run into potential bugs)")
    parser.add_argument('-sample_id','--sample_id',dest="sample_id",help="CHGVID",required=True)
    parser.add_argument('-block_size','--block_size',dest="B",type=int,default=10000,help="The Block Size to be used for binning(default 10000)")
    parser.add_argument('-output_dir','--output_dir',dest="output_dir",help="Directory to save the output files(will be created if it does not exist)",required=True)
    parser.add_argument('inputfile',metavar='bed/gvcf file',help="The file to bin(can be gzipped)")
    parser.add_argument('-mode','--mode',choices=['gq','coverage'])
    args = parser.parse_args()
    
    main(args)

