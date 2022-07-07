# dragen_pipe

To run the gatk pipeline : 

change directory to dragen/python 

`luigi --module gatk_pipe ArchiveSample --pseudo-prepid {pseudo_prepid} --workers {workers}`<br> where pseudo_prepid is the pseudo_prepid for a sample and workers is the number of workers to use (normally 2).

To run automation script:<br> 
`./run_gatk.py`<br>
--help will list standard parameters, e.g. to adjust the number of samples to process concurrently, and the number of workers per job.  --additional_sample_requirements may be used to customize the query used to get samples (it simply appends anything added here to the query, see dragen_db_statements.py.GET_SAMPLES for the query).  This script will also accept any arbitrary parameters not listed here and pass them on unchanged to luigi, so look at code/class definitions/luigi.cfg as appropriate/needed.

To run manually using parallel :

cat {sample_file} | parallel -j {n_jobs} --joblog {parallel_log_file} --tmpdir /nfs/seqscratch_ssd/rp2801/parallel_tmp/ /nfs/goldstein/software/python2.7.7/bin/python run_gatk_jobs.py

To import samples:<br>
ssh to annodb02<br>
cd to import/src<br>
`./import_samples.py -d waldb4 --run_locally`<br>
Additional parameters may be viewed via --help.
