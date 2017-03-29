#!/usr/bin/awk -f


BEGIN{
    FS="\t"
}

{
    ## Need to have a convoluted logic here because some of the production vcf headers have white spaces instead
    ## of tabs in the sample header line only , here is an e.g. :
    ## /nfs/seqsata10/ALIGNMENT/BUILD37/EXOME/GATK/isnd24071k2/combined/isnd24071k2.analysisReady.annotated.vcf.gz
    ## look for the line starting with #CHROM by splitting it by a tab in python or your favorite language
    
    if($0 ~/\#/) {
	if ($0 ~/\#CHROM/){
	    ## the variable $sample is passed from the command line and refers to the sample name
	    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample";	    
	}
	else
	    print
    }
    else if ($7 == "PASS")
	print
}
END{}
