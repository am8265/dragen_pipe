vcf=$1
tempfile=$2

cat ${vcf} | awk 'NR==9{print "##INFO=<ID=SAO,Number=1,Type=Integer,Description=\"Variant Allele Origin: 0. unspecified, 1. Germline, 2. Somatic, 3. Both\">"}1' > ${tempfile}

cmd="mv ${tempfile} ${vcf}"
echo $cmd
eval $cmd

#cmd="bgzip -f ${vcf}"
#echo $cmd
#eval $cmd

#cmd="/nfs/goldstein/software/bin/tabix -f ${vcf}.gz"
#echo $cmd
#eval $cmd
