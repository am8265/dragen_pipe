# dragen_pipe

To run the gatk pipeline : 

change directory to dragen/python 

luigi --module gatk_pipe ArchiveSample --sample-name {the chgvid} --pseudo-prepid {the pseudo_prepid} --capture-kit-bed {the capture bed file} --sample-type {exome/genome/custom capture} --scratch {the scratch directory} --dont-remove-tmp-dir --poll-time 120 --workers 2 --worker-wait-interval 120 --scheduler-remove-delay 86400

To run automation script : 
python gatk_wrapper.py --max-processes 300 --wait-time 1000 

To run manually using parallel :

cat {sample_file} | parallel -j {n_jobs} --joblog {parallel_log_file} --tmpdir /nfs/seqscratch09/rp2801/parallel_tmp/ /nfs/goldstein/software/python2.7.7/bin/python run_gatk_jobs.py
