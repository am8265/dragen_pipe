
import os

for i in range(1):
    print("WORKING ON CASE {}".format(i)) 
    cmd = "/nfs/goldstein/software/python2.7.7/bin/luigi  --module ethnicity_check_dragen PredictAndUpdate --base-output-directory  /nfs/seqscratch_ssd/tmp --pseudo-prepid 0 --local-scheduler"
    os.system(cmd)
 
