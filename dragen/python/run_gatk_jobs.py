import sys
import os
from datetime import datetime

info=sys.argv[1]

pseudo_prepid,chgvid,seqtype,capturekit_bed,scratch = info.split('\t')
now=datetime.now()
timestamp=now.strftime("%Y%m%d")


out_dir=os.path.join("/nfs/seqscratch09/rp2801/luigi_logs/",timestamp)

## Check if the directory exists : 
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)
error_file=os.path.join(out_dir,chgvid+'.'+pseudo_prepid+'.stderr')
out_file=os.path.join(out_dir,chgvid+'.'+pseudo_prepid+'.stdout')

cmd="""luigi --module gatk_pipe ArchiveSample --sample-name {0} --pseudo-prepid {1} --capture-kit-bed  {2} --sample-type {3} --scratch {4} --poll-time 30 --workers 2 --worker-wait-interval 120 --scheduler-remove-delay 86400 > {5} 2>{6}""".format(chgvid,pseudo_prepid,capturekit_bed,seqtype,scratch,out_file,error_file)
#cmd="""luigi --module gatk_pipe GenotypeGVCFs --sample-name {0} --pseudo-prepid {1} --capture-kit-bed  {2} --sample-type {3} --dont-remove-tmp-dir --no-tarball --poll-time 120 --workers 2 --worker-wait-interval 180 --scheduler-remove-delay 86400 > {4} 2>{5}""".format(chgvid,pseudo_prepid,capturekit_bed,seqtype,out_file,error_file)

os.system(cmd)
    
