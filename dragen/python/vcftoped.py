import os,sys,subprocess,string
from optparse import OptionParser
import argparse
import os
import gzip

# read sample id, sample fid, vcf

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

parser = OptionParser()
parser.add_option("--scratch", type = "string", dest = "scratch", default = ".")
parser.add_option("--vcffile", type = "string", dest = "vcffile")
parser.add_option("--markers", type = "string", dest = "markerfile")
parser.add_option("--refs", type = "string", dest = "reffile")
parser.add_option("--sample", type = "string",dest = "sample")
(options,args) = parser.parse_args()


markerfile = open(options.markerfile).readlines()#map file
markers = []
for line in markerfile:
    temp = line.split("\t")
    markers.append(temp[0]+"."+temp[-1].strip("\n"))

reffile = open(options.reffile).readlines()
refs = []
geno = []
for line in reffile:
    refs.append(line[0]+line[1])
    geno.append("11")
    

#vcffile = open(options.vcffile).readlines()
vcffile = fh(options.vcffile)

for line in vcffile:
    if line[0] != "#":
        temp = line.split("\t")
        chrom = temp[0]
        pos = temp[1]
        
        if markers.count(chrom+"."+pos) == 1:
            site = markers.index(chrom+"."+pos)
            allele1 = temp[3]
            allele2 = temp[4]
            if refs[site].count(allele1) == 1 and refs[site].count(allele2) == 1:# check triallic
                depth1 = line.split("DP=")[1]
                depth = string.atoi(depth1.split(";")[0])
               
                if depth >= 10:
                    AN = line.split("AN=")[1][0]
                    AC = line.split("AC=")[1][0]
                    if AN == "2" and AC == "0":
                        geno[site] = allele1 + allele1
                    elif AN == "2" and AC == "1":
                        geno[site] = allele1 + allele2
                    elif AN == "2" and AC == "2":
                        geno[site] = allele2 + allele2
                else:
                    geno[site] = "00"

            else:
                geno[site] = "00"
                
    if line[0] == "X":
        break

outfile = open(options.scratch+options.sample+".geno1","w")
for line in geno:
    outfile.write(line+"\n")

outfile.close()
