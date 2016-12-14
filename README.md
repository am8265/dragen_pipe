# dragen_pipe

To run the gatk pipeline : 

change directory to dragen/python 

luigi --module gatk_pipe ArchiveSample --sample-name {the chgvid} --pseudo-prepid {the pseudo_prepid} --capture-kit-bed {the capture bed file} --sample-type {exome/genome/custom capture} --scratch {the scratch directory} --dont-remove-tmp-dir --poll-time 120 --workers 2 --worker-wait-interval 120 --scheduler-remove-delay 86400
