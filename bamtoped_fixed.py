import os,sys,string,subprocess
from optparse import OptionParser

def getcoverage(site,bamfile,scratch,sample,cluster):
    chromosome = site.split(".")[0]
    site = string.atoi(site.split(".")[1])
    if (cluster == "lsrc" or cluster == "LSRC"):
        samtool = "/nfs/goldstein/goldsteinlab/software/samtools-0.1.17/samtools view "+bamfile+" "+chromosome+":"+str(site)+"-"+str(site+1)+ " -bh > "+scratch+sample+".temp.bam"
    else:
        samtool = "/nfs/chgv/seqpipe01/SOFTWARE/samtools-0.1.17/samtools view "+bamfile+" "+chromosome+":"+str(site)+"-"+str(site+1)+ " -bh > "+scratch+sample+".temp.bam"
        
    subprocess.call(samtool,shell=True)

    FNULL = open(os.devnull, 'w')
    #print bamfile
    if (cluster == "lsrc" or cluster == "LSRC"):
        subprocess.call("/nfs/goldstein/goldsteinlab/software/samtools-0.1.17/samtools mpileup "+scratch+sample+".temp.bam > "+scratch+sample+".temp.pileup",shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    else:
        subprocess.call("/nfs/chgv/seqpipe01/SOFTWARE/samtools-0.1.17/samtools mpileup "+scratch+sample+".temp.bam > "+scratch+sample+".temp.pileup",shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        
    lines = open(scratch+sample+".temp.pileup").readlines()
    depth = "0"
    refbase = "0"
    
    for line in lines:
        temp = line.split("\t")
        if temp[0] == chromosome and temp[1] == str(site):
            depth = temp[3]
            bases = temp[4].lower()
            numA = bases.count("a")
            numC = bases.count("c")
            numG = bases.count("g")
            numT = bases.count("t")

            if (numA > numG and numA > numC and numA > numT):
                refbase = "A"
            if (numC > numG and numC > numA and numC > numT):
                refbase = "C"
            if (numG > numA and numG > numC and numG > numT):
                refbase = "G"
            if (numT > numG and numT > numC and numT > numA):
                refbase = "T"
        
    subprocess.call("rm "+scratch+sample+".temp.pileup "+scratch+sample+".temp.bam",shell=True)

    return string.atoi(depth),refbase


parser = OptionParser()
parser.add_option("--scratch", type = "string", dest = "scratch", default = ".")
parser.add_option("--bamfile", type = "string", dest = "bamfile")
parser.add_option("--markers", type = "string", dest = "markerfile")
parser.add_option("--refs", type = "string", dest = "reffile")
parser.add_option("--sample", type = "string",dest = "sample")
parser.add_option("--cluster",type="string",dest = "cluster",default="lsrc")
(options,args) = parser.parse_args()

if options.scratch[-1] != "/":
    options.scratch = options.scratch + "/"

markerfile = open(options.markerfile).readlines()
markers = []
for line in markerfile:
    temp = line.split("\t")
    markers.append(temp[0]+"."+temp[-1].strip("\n"))

reffile = open(options.reffile).readlines()
genofile = open(options.scratch + options.sample+".geno1").readlines()
refs = []
geno = []
for line in reffile:
    refs.append(line.strip("\n"))
for line in genofile:
    geno.append(line.strip("\n"))


for site in range(len(geno)):
    if geno[site] == "11": #none at vcf
        # bam file check
        depth,refbase = getcoverage(markers[site],options.bamfile,options.scratch,options.sample,options.cluster)
        if depth >= 10:
            geno[site] = refbase + refbase

                      
        else:
            geno[site] = "00"

outfile = open(options.scratch + options.sample+".geno.txt","w")

for line in geno:
    outfile.write(line + "\n")
outfile.close() 
