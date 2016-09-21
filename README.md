# dragen_pipe

# To run the gatk pipeline : 
# change directory to dragen/python 
# luigi --module gatk_pipe ArchiveSample --sample-name <the chgvid> --pseudo-prepid <the pseudo_prepid> --capture-kit-bed <the capture bed file> --sample-type <exome/genome/custom capture> --dont-remove-tmp-dir --no-tarball --poll-time 120 --workers 3 --worker-wait-interval 180 --scheduler-remove-delay 86400
