import sys
import subprocess
import argparse
import os


#######################################
###                                 ###
###  Script for binning GQ values   ###
###  from a gvcf file               ###
###  Author : Raghav                ###
#######################################

def is_gzipped(file_name):
    """is the specified file gzipped?
    """
    with open(file_name, "rb") as fh:
        magic_number = fh.read(2)
    return magic_number == "\x1f\x8b"

def fh(file_name):
    """return a file handle to the file whether it's gzipped or not
    """
    if os.path.isfile(file_name):
        if is_gzipped(file_name):
            return gzip.open(file_name)
        else:
            return open(file_name)
    else:
        raise argparse.ArgumentTypeError(
            "{} does not exist".format(file_name))

def gq_to_letter(gq):
    
    """ Returns Letter encoding for gq value
        gq => numerical gq value
    """

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

    
def bin_gq_to_letter(gq_val, gq_len):
    """ Return the printable form for the alpha-numeric gq bins
        i.e. gq_val = ['a','b','c'] gq_len = [200,300,500]
        200a,300b,500c is returned
        gq_val => a list of letter encoded gq values
        gq_len => a list of integers corresponding to the lengths of the
                        above letter encoded gq values
    """

    return_gq_string_list = []
    check_len = 0
    for c,l in zip(gq_val, gq_len):
        check_len+=l
        l = radix36_encoding(l)
        return_gq_string_list.append(str(l)+str(c))
        
    if check_len != 10000:
        raise Exception("Incorrect Binning String Length!")
    
    return ''.join(return_gq_string_list)


def update_gq_bins(n, gq, gq_len, gq_val):
    
    """ Update the coded gq values according to the letter encodings
        i.e. if n = 100 , gq = 'a' , gq_val = ['b','c','a'] and 
        gq_len = [100,200,200], the program will update the gq_len
        to [100,200,300] and gq_val to ['b','c','a']
        n => numerical value of gq
        gq => letter code for gq bin
        gq_len => list of numerical gq values seen in the current block
        gq_val => list of letter codes for gq bins seen in the current block
    """

    if(not gq_len):  # If list is empty
        gq_len.append(n)
        gq_val.append(gq)

    else:
        if(gq == gq_val[-1]):
            gq_len[-1] = gq_len[-1] + n
            
        else:
            gq_len.append(n)
            gq_val.append(gq)

    return [gq_len, gq_val]

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
        raise Exception("Negative bin length encountered !")
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
    line = line.strip('\n')

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
            print line 
        end = start

    

    return [chromosome,start,end,gq]


def is_empty_gqlist(gq_len):
    """
    Check whether the gq_bin length list is empty
    returns a boolean
    """
    if gq_len:
        return False 
    else:
        return True

def output_line(gq_str,exclusion_str,sample_id,out_start,block,bin_out):
    """
    Output gq_str to a file as long as it is not the
    exclusion string
    """

    if gq_str != exclusion_str:
        bin_out.write("%s\t%s\n" % (
                str(out_start/block),gq_str))
    
        
def chromosome_check(prev_chromosome,chromosome):
    """
    Check whether the chromsome has changed 
    """

    return prev_chromosome == chromosome 
                
    
    
def reinitialize_variables(out_start,block):
    """
    Reinitializes key variables after the a block has been spanned
    """

    gq_len = []
    gq_val = []
    out_start = out_start + block
    current_block = block
    
    return (out_start,current_block,gq_len,gq_val)


def main():
    """ Main Function """

    block = int(sys.argv[1])
    gvcf = sys.argv[2]
    sample_id = sys.argv[3]
    out_dir = sys.argv[4]
    
    IN = fh(gvcf)
    
    
    human_chromosomes = range(1,23)
    human_chromosomes = [str(x) for x in human_chromosomes]
    human_chromosomes.extend(['X','Y','MT'])
    prev_chromosome = '1'
    exclusion_str = radix36_encoding(block) + 'a'
    flag = 0

    gq_len = []
    gq_val = []
    run_len = []
    
    for line in IN:
        line = line.strip('\n')
        if line[0] == '#': ## skip header
            continue
        
        chromosome,start,end,gq = parse_gvcf_gq(line)
        
        

        if chromosome not in human_chromosomes:
            continue ## Skip scaffolds
        
        if chromosome_check(prev_chromosome,chromosome) == False:
            if is_empty_gqlist(gq_len) == False: ## Account for remaining bins if present
                if block - sum(gq_len) != 0:  ## Make sure the block is filled
                    update_gq_bins(
                        block - sum(gq_len),'a',gq_len,gq_val)
                
                    gq_str = bin_gq_to_letter(
                        gq_val,gq_len)
    
                    output_line(gq_str,exclusion_str,sample_id,out_start,block,bin_out)
                    out_start,current_block,gq_len,gq_val = reinitialize_variables(out_start,block)
            bin_out.close()
            #bin_out = open(os.path.join(out_dir,sample_id+'_gq_binned_'+str(block)+'_chr%s.txt'%chromosome),'w')
            flag = 0 
                

        if flag == 0:
            bin_out = open(os.path.join(out_dir,sample_id+'_gq_binned_'+str(block)+'_chr%s.txt'%chromosome),'w')
            out_start = start - start % block
            gap = start % block - 1
            current_block = block
            gq_len = []
            gq_val = []
        else:
            gap = start - prev_stop - 1 

        if (gap > 0): # Take account of gaps, if they are present
            if(current_block - gap > 0): ## block spans the gap
                current_block = current_block - gap
                n = gap
                gq_letter = 'a'
                gq_len,gq_val = update_gq_bins(
                    n,gq_letter,gq_len,gq_val)
                
            elif(current_block - gap == 0): ## The gap exactly spans the block
                n = gap
                gq_letter = 'a'
                gq_len,gq_val = update_gq_bins(
                    n,gq_letter,gq_len,gq_val)
                ## Get printable form for coded gq
                gq_str = bin_gq_to_letter(
                gq_val,gq_len)
                output_line(gq_str,exclusion_str,sample_id,out_start,block,bin_out)
                ## Reinitialization
                out_start,current_block,gq_len,gq_val = reinitialize_variables(
                out_start,block)
                
            elif (current_block - gap < 0): ## The gap spans the block
                n = current_block
                gq_letter = 'a'
                gq_len,gq_val = update_gq_bins(
                    n,gq_letter,gq_len,gq_val)
                gq_str = bin_gq_to_letter(
                    gq_val,gq_len)
                output_line(gq_str,exclusion_str,sample_id,out_start,block,bin_out)
                ## Update remaining gap
                remaining_gap = gap - current_block
                ## Reinitialization
                out_start,current_block,gq_len,gq_val = reinitialize_variables(
                    out_start,block)

                num_blocks = remaining_gap/block
                # Check to see if the gap still spans another block
                if(current_block - remaining_gap <= 0):
                    ## Change remaining gap to be the residual gap
                    remaining_gap = remaining_gap % current_block
                        
                    while(num_blocks > 0): # Iterate over the blocks
                        n = current_B
                        gq_letter = 'a'
                        gq_len,gq_val = update_gq_bins(
                            n,gq_code,gq_len,gq_val)
                        gq_str = bin_gq_to_letter(
                            gq_val,gq_len)
                        output_line(gq_str,exclusion_str,sample_id,out_start,block,bin_out)
                        num_blocks = num_blocks - 1
                        ## Reinitialization
                        out_start,current_block,gq_len,gq_val = reinitialize_variables(
                            out_start,block)
                        
                    if(remaining_gap > 0):
                        n = remaining_gqp
                        gq_letter = 'a'
                        gq_len,gq_val = update_gq_bins(
                            n,gq_letter,gq_len,gq_val)
                        # Update block size
                        current_block = current_block - remaining_gap
                else: ## If the remaining gap is lesser than a block size
                    n = remaining_gap
                    gq_letter = 'a'
                    gq_len,gq_val = update_gq_bins(
                        n,gq_letter,gq_len,gq_val)
                    # Update block size
                    current_block = current_block - remaining_gap

        ## Accounting for intervals
        current_interval = int(end) - int(start) + 1      
        current_gq = int(gq)
        gq_letter = gq_to_letter(current_gq)

        if current_block - current_interval > 0:
            # Update block size
            current_block = current_block - current_interval
            n = current_interval
            gq_len,gq_val = update_gq_bins(
                n,gq_letter,gq_len,gq_val)

        elif current_block - current_interval == 0:
            n = current_interval
            gq_len,gq_val = update_gq_bins(
                n,gq_letter,gq_len,gq_val)
            gq_str = bin_gq_to_letter(
                gq_val,gq_len)

            ## print output to file
            
            output_line(gq_str,exclusion_str,sample_id,out_start,block,bin_out)
        
            ## Reinitialization
            out_start,current_block,gq_len,gq_val = reinitialize_variables(
                out_start,block)
           
                
        elif current_block - current_interval < 0:
            n = current_block
            gq_len,gq_val = update_gq_bins(
                n,gq_letter,gq_len,gq_val)

            gq_str = bin_gq_to_letter(
                gq_val,gq_len)
            
            output_line(gq_str,exclusion_str,sample_id,out_start,block,
                        bin_out)

            ## Reinitialization
            out_start,current_block,gq_len,gq_val = reinitialize_variables(
                out_start,block)            
            remaining_interval = current_interval - n
            num_blocks = remaining_interval/block
            
            
            if(current_block - remaining_interval <= 0):
                remaining_interval = remaining_interval % current_block
                while(num_blocks > 0):
                    n = current_block
                    gq_len, gq_val = update_gq_bins(
                        n, gq_letter, gq_len, gq_val)
                    gq_str = bin_gq_to_letter(
                        gq_val, gq_len)
                    
                    output_line(gq_str,exclusion_str,sample_id,out_start,block,
                        bin_out)
                    
                    ## Reinitialization
                    out_start,current_block,gq_len,gq_val = reinitialize_variables(
                        out_start,block)
                    num_blocks = num_blocks - 1
                    

                if(remaining_interval > 0):
                    n = remaining_interval
                    gq_len, gq_val = update_gq_bins(
                        n, gq_letter, gq_len, gq_val)
                    # Update block size
                    current_block = current_block - remaining_interval

            else:
                n = remaining_interval
                gq_len, gq_val = update_gq_bins(
                    n, gq_letter, gq_len, gq_val)
                # Update block size
                current_block = current_block - remaining_interval
                  
        flag = 1
        # Store previous values
        prev_start = start
        prev_stop = end
        prev_chromosome = chromosome
        
        
    if gq_len : ## Account for remaining bins if present
        # Add a's if the block is not fully spanned
        sum_length = sum(gq_len)
        if block - sum_length != 0:
            gq_len,gq_val = update_gq_bins(block - sum_length, 'a',
                                           gq_len,gq_val)


        gq_str = bin_gq_to_letter(gq_val,gq_len)
        output_line(gq_str,exclusion_str,sample_id,out_start,block,
                        bin_out)
            
            
    IN.close()
    bin_out.close()

    
if __name__ == "__main__":
    main()
