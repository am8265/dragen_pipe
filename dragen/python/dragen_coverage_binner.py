import sys
import os
import subprocess
from utilities import is_gzipped,fh

###############################################################################
###       Script to convert dragen bed file to binned coverage format       ###
###       Run as : "Run as : python/pypy dragen_coverage_binner.py          ###
###       <Block_Size> <Sample_ID> <Bed_File> <Output_Dir>                  ###
###                                                                         ###
###       Author : Raghav                                                   ###
###       Note : For faster run times, use pypy instead of standard python  ###
###              interpreter                                                ###
###############################################################################


def get_sample_id(dragen_file):
    
    """ Returns the sample id to print in the output file
    """
    
    return dragen_file.split('/')[-1].split('.')[0]


def coverage_to_letter(coverage):
    
    """ Returns Letter encoding for coverage value
        coverage => numerical coverage value
    """

    if(coverage < 3):
        return 'a'
    elif(coverage < 10):
        return 'b'
    elif(coverage < 20):
        return 'c'
    elif(coverage < 30):
        return 'd'
    elif(coverage < 40):
        return 'e'
    elif(coverage < 50):
        return 'f'
    elif(coverage < 60):
        return 'g'
    elif(coverage < 80):
        return 'h'
    elif(coverage < 100):
        return 'i'
    elif(coverage < 150):
        return 'j'
    elif(coverage < 200):
        return 'k'
    elif(coverage < 250):
        return 'l'
    elif(coverage < 300):
        return 'm'
    elif(coverage < 500):
        return 'n'
    elif(coverage < 600):
        return 'o'
    elif(coverage < 800):
        return 'p'
    elif(coverage < 1000):
        return 'q'
    else:
        return 'r'


def bin_coverage_to_letter(coverage_val, coverage_len):
    
    """ Return the printable form for the alpha-numeric coverage bins
        i.e. coverage_val = ['a','b','c'] coverage_len = [200,300,500]
        200a,300b,500c is returned
        coverage_val => a list of letter encoded coverage values
        coverage_len => a list of integers corresponding to the lengths of the
                        above letter encoded coverage values
    """

    return_coverage_string_list = []
    check_len = 0
    for c,l in zip(coverage_val, coverage_len):
        check_len+=l
        l = radix36_encoding(l)
        return_coverage_string_list.append(l+c) 

    if check_len!=10000:
        raise Exception("Incorrect binning string length")
    
    return ''.join(return_coverage_string_list)


def update_coverage_bins(n, coverage, coverage_len, coverage_val):
    
    """ Update the coded coverage values according to the letter encodings
        i.e. if n = 100 , coverage = 'a' , coverage_val = ['b','c','a'] and 
        coverage_len = [100,200,200], the program will update the coverage_len
        to [100,200,300] and coverage_val to ['b','c','a']
        n => numerical value of coverage
        coverage => letter code for coverage bin
        coverage_len => list of numerical coverage values seen in the current block
        coverage_val => list of letter codes for coverage bins seen in the current block
    """

    if(not coverage_len):  # If list is empty
        coverage_len.append(n)
        coverage_val.append(coverage)

    else:
        if(coverage == coverage_val[-1]):
            coverage_len[-1] = coverage_len[-1] + n

        else:
            coverage_len.append(n)
            coverage_val.append(coverage)

    return [coverage_len, coverage_val]

def radix36_encoding(number):
    
    """ 
    number => numerical coverage value in decimal
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
        raise Exception("Negative bin length encoutered, bug in your code !")

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


def main():
    """ Main Function """

    if len(sys.argv) != 5:
        print "Invalid Command Line Arguments !"
        print ("Run as : python/pypy dragen_coverage_binner.py <Block_Size>"
               "<Sample_ID> <Bed_File> <Output_Dir>")
        sys.exit(1)

    else:
        # Command Line Arguments
        B = int(sys.argv[1])
        sample_id = sys.argv[2]
        dragen_bed_file = sys.argv[3]
        output_dir = sys.argv[4]
        
        # Initialize variables and lists
        flag = 0
        coverage_len = []
        coverage_val = []
        human_chromosomes = []
        current_B = B
        #sample_id = get_sample_id(dragen_bed_file)
        prev_chromosome = '1'
        exclusion_str = radix36_encoding(B) + 'a'

        # Store Human Chromosome Numbers and Letters, will be used to
        # filter out scaffolds
        human_chromosomes.extend(range(1, 23))
        human_chromosomes = [str(x) for x in human_chromosomes]
        human_chromosomes.extend(['X', 'Y','MT'])

        # Open File Handles For Input and Output Files
        dragen_bed_in = fh(dragen_bed_file)
        
        bin_out = open(os.path.join(output_dir,sample_id + "_coverage_binned_" + str(B) + "_chr1.txt"), 'w')

        for line in dragen_bed_in:
            line = line.strip('\n')
            contents = line.split('\t')
            chromosome = contents[0]
            start = int(contents[1])
            stop = int(contents[2])
            coverage = int(contents[3])

            ######## INITIALIZATION SECTION #######
            
            # Convert 0-based to 1-based
            start = start + 1
            stop = stop + 1

            if(chromosome not in human_chromosomes):
                continue      # Skip iterations for scaffolds

            if(prev_chromosome != chromosome):
                if(coverage_len):  # Account for remaining bins if present
                    # Add a's if the block is not fully spanned
                    sum_length = sum(coverage_len)
                    if(B - sum_length != 0):
                        update_coverage_bins(
                            B - sum_length, 'a', coverage_len, coverage_val)

                    coverage_str = bin_coverage_to_letter(
                        coverage_val, coverage_len)
                    if(coverage_str != exclusion_str):
                        bin_out.write("%s\t%s\n" % (
                            str(out_start/B), coverage_str))

                # File handles are changed when chromosome changes
                bin_out.close()
                bin_out = open(os.path.join(output_dir,sample_id + "_coverage_binned_" + str(B) + "_chr%s.txt" % chromosome), 'w')
                flag = 0

            if(flag == 0):  # Only for the first line of a new chromosome
                out_start = start - start % B + 1
                gap = start % B - 1
                current_B = B
                coverage_len = []
                coverage_val = []
            else:
                gap = start - prev_stop
                #print gap,out_start,coverage,current_B
                
            ######### END OF INITALIZATION SECTION ##########

            ################ ACTUAL LOGIC PROCESSING STARTS HERE ############
            ## Divided into two stages, first is to process any gaps and then
            ## proceed with processing the actual interval
            
            ## gaps occur in the event that the start position begins at an arbitrary
            ## position skipping any genome coordinates from the beginning of a block
            ## all of these bases are assinged the bin coverage 'a'

            ###### FIRST process the gaps ######
            if(gap > 0):  # Take account of gaps, if they are present
                if(current_B - gap > 0):  # Current block size remaining spans
                                          # the gap
                    current_B = current_B - gap  # Account for reduction in
                    # current block size
                    n = gap
                    coverage_code = 'a'  # Note : gaps are coded as 'a's
                    coverage_len, coverage_val = update_coverage_bins(
                        n, coverage_code, coverage_len, coverage_val)

                elif(current_B - gap == 0):  # The gap exactly spans the block
                                             # size
                    n = gap  # Numerical Length of the coded coverage
                    coverage_code = 'a'
                    coverage_len, coverage_val = update_coverage_bins(
                        n, coverage_code, coverage_len, coverage_val)  # Update bins
                    # Get the printable output for coded coverage
                    coverage_str = bin_coverage_to_letter(
                        coverage_val, coverage_len)
                    if(coverage_str != exclusion_str):
                        bin_out.write("%s\t%s\n" % (
                            str(out_start/B), coverage_str))

                    # Reinitialize key variables
                    current_B = B
                    gap = 0
                    coverage_len = []
                    coverage_val = []
                    out_start = out_start + B

                elif(current_B - gap < 0):  # The gap spans the current block
                                            # size
                    n = current_B
                    coverage_code = 'a'
                    coverage_len, coverage_val = update_coverage_bins(
                        n, coverage_code, coverage_len, coverage_val)
                    coverage_str = bin_coverage_to_letter(
                        coverage_val, coverage_len)
                    if(coverage_str != exclusion_str):
                        bin_out.write("%s\t%s\n" % (
                            str(out_start/B), coverage_str))

                    # Update remaining gap
                    remaining_gap = gap - current_B
                    
                    # Reinitialize key variables
                    current_B = B
                    out_start = out_start + B
                    num_blocks = remaining_gap / B
                    coverage_len = []
                    coverage_val = []

                    # Checking to see if the gap still spans another block size
                    if(current_B - remaining_gap <= 0):
                        # Change remaining gap to be the residual gap left
                        # after dividing the gaps into discrete block steps
                        remaining_gap = remaining_gap % current_B
                        
                        
                        while(num_blocks > 0):  # Iterate over the number of
                                                # blocks
              
                            n = current_B
                            coverage_code = 'a'
                            coverage_len, coverage_val = update_coverage_bins(
                                n, coverage_code, coverage_len, coverage_val)
                            coverage_str = bin_coverage_to_letter(
                                coverage_val, coverage_len)
                            if(coverage_str != exclusion_str):
                                bin_out.write("%s\t%s\n" % (
                                    str(out_start/B), coverage_str))
                            num_blocks = num_blocks - 1

                            # Reinitialize key variables
                            coverage_len = []
                            coverage_val = []
                            out_start = out_start + B
                            current_B = B
                            
                        if(remaining_gap > 0):
                            n = remaining_gap
                            coverage_code = 'a'
                            coverage_len, coverage_val = update_coverage_bins(
                                n, coverage_code, coverage_len, coverage_val)
                            # Update block size
                            current_B = current_B - remaining_gap

                    else:  # If the remaining gap is lesser than a block size
                        n = remaining_gap
                        coverage_code = 'a'
                        coverage_len, coverage_val = update_coverage_bins(
                            n, coverage_code, coverage_len, coverage_val)

                        # Update block size                
                        current_B = current_B - remaining_gap
            ########## Finished processing gaps #######
                        
            ###### Begin Processing the actual interval now ######
            interval_len = stop - start
            if(current_B - interval_len > 0):  # Current block size
                                               # spans the interval
                # Update block size
                current_B = current_B - interval_len
                n = interval_len
                coverage_code = coverage_to_letter(coverage)
                coverage_len, coverage_val = update_coverage_bins(
                    n, coverage_code, coverage_len, coverage_val)

            elif(current_B - interval_len == 0):  # Current block size exactly
                                                  #  spans the interval
                n = interval_len
                coverage_code = coverage_to_letter(coverage)
                coverage_len, coverage_val = update_coverage_bins(
                    n, coverage_code, coverage_len, coverage_val)
                coverage_str = bin_coverage_to_letter(
                    coverage_val, coverage_len)
                if(coverage_str != exclusion_str):
                    bin_out.write("%s\t%s\n" %
                                  (str(out_start/B), coverage_str))

                # Reinitialize key variables
                out_start = out_start + B
                current_B = B
                coverage_len = []
                coverage_val = []

            elif(current_B - interval_len < 0):  # The interval spans current
                                                 # block size
                n = current_B
                coverage_code = coverage_to_letter(coverage)
                coverage_len, coverage_val = update_coverage_bins(
                    n, coverage_code, coverage_len, coverage_val)
                coverage_str = bin_coverage_to_letter(
                    coverage_val, coverage_len)
                if(coverage_str != exclusion_str):
                    bin_out.write("%s\t%s\n" %
                                  (str(out_start/B), coverage_str))
                out_start = out_start + B

                #remaining_interval = interval_len - current_B
                remaining_interval = interval_len - n
                # Reinitialize key variables
                current_B = B
                num_blocks = remaining_interval / B
                coverage_len = []
                coverage_val = []

                # Similar logic as was used for the gaps before
                # Iterate over discrete blocks and do the coverage
                # updates and any left over residual interval 
                if(current_B - remaining_interval <= 0):
                    remaining_interval = remaining_interval % current_B
                    while(num_blocks > 0):
                        n = current_B
                        coverage_code = coverage_to_letter(coverage)
                        coverage_len, coverage_val = update_coverage_bins(
                            n, coverage_code, coverage_len, coverage_val)
                        coverage_str = bin_coverage_to_letter(
                            coverage_val, coverage_len)
                        if(coverage_str != exclusion_str):
                            bin_out.write("%s\t%s\n" % (
                                str(out_start/B), coverage_str))
                        num_blocks = num_blocks - 1

                        # Reinitialize key variables
                        coverage_len = []
                        coverage_val = []
                        out_start = out_start + B
                        current_B = B

                    if(remaining_interval > 0):
                        n = remaining_interval
                        coverage_code = coverage_to_letter(coverage)
                        coverage_len, coverage_val = update_coverage_bins(
                            n, coverage_code, coverage_len, coverage_val)
                        # Update block size
                        current_B = current_B - remaining_interval

                else:
                    n = remaining_interval
                    coverage_code = coverage_to_letter(coverage)
                    coverage_len, coverage_val = update_coverage_bins(
                        n, coverage_code, coverage_len, coverage_val)
                    # Update block size
                    current_B = current_B - remaining_interval

            flag = 1
            # Store previous values
            prev_start = start
            prev_stop = stop
            prev_chromosome = chromosome

        if(coverage_len):  # Account for remaining bins if present
            # Add a's if the block is not fully spanned
            sum_length = sum(coverage_len)
            if(B - sum_length != 0):
                coverage_len,coverage_val = update_coverage_bins(B - sum_length, 'a',
                                     coverage_len, coverage_val)

            coverage_str = bin_coverage_to_letter(coverage_val, coverage_len)
            if(coverage_str != exclusion_str):
                bin_out.write("%s\t%s\n" %
                              (str(out_start/B), coverage_str))

        dragen_bed_in.close()
        bin_out.close()


if __name__ == "__main__":
    main()
