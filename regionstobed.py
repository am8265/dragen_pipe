import sys

def output_bed(regions_file,output_bed):
    """ Parse the atav regions bed format
    and output to a regular bed file 
    """

    with open(regions_file,'r') as IN, open(output_bed,'w') as OUT:
        for line in IN:
            contents = line.strip('\n').split(' ') 
            gene = contents[0]
            chromosome = contents[1]
            col3 = contents[2]
            ## Remove the terminal brackets 
            col3 = col3.strip('(')
            col3 = col3.strip(')')
            ## Process the different regions which were inside the brackets that are seperated by a ',' 
            coordinates_col3 = col3.split(',')
            for coordinates in coordinates_col3: ## Iterate over the exonic coordinates ? 
                start,stop = coordinates.split('..') ## The start,stop are obtained
                OUT.write("%s\t%s\t%s\t%s\n"%(chromosome,start,stop,gene))






