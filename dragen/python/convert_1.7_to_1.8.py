#!/nfs/goldstein/software/python2.7.7/bin/python
"""Converts Ilumina 1.7 fastq sequencing performed at Duke University via CHGV
into Illumina 1.8 fastq format.
"""

import shlex
import subprocess
import sys
import tarfile
from glob import glob

def main():
    sample_name = sys.argv[1]
    sample_type = sys.argv[2]
    global run_script
    run_script = 0
    counter = 1

    for tar_file in glob('/nfs/fastq16/{}_1.7/{}/*tar.gz'.format(sample_type.upper(),sample_name)):
        dest = '/nfs/qumulo/{}/{}/{}/'.format(sample_type.upper(),sample_name,counter)
        read_3_flag = 0

        mkdir_cmd = 'mkdir -p {}'.format(dest)
        runCmd(mkdir_cmd)

        if run_script:
            print "Extracting {}".format(tar_file)
            tfile = tarfile.open(name=tar_file, mode = 'r:gz')
            tfile.extractall(path = dest)
            tfile.close()

        sequence_files = glob('{}/s*_sequence.txt'.format(dest))
        if len(sequence_files) == 0:
            raise Exception, "No sequencing files were found"

        for sequence_file in sequence_files:
            info = sequence_file.split('_')
            """All Casava 1.7 samples should have their fastq files named in
            the following schema: s_[lane]_[read_number]_sequence.txt.

            Ex: s_2_1_sequence.txt

            Due to this naming schema there should be no more than two files
            (one for each read) on a lane so there is no need to increment the
            '001' string at the end of the new file name.

            In some rare cases read 2 was numbered as read 3 following
            Illumina's method of naming the reads.  Added a check for this
            so that "read 2" or the index read is not present and converted
            """
            if info[-1] != 'sequence.txt':
                raise Exception, 'File {} is not in the proper naming format!'.format(sequence_file)

            if info[-2] == '1':
                new_file_name = '{}_AAAAAA_L00{}_R{}_001.fastq.gz'.format(sample_name,info[-3],'1')
            elif info[-2] == '3':
                new_file_name = '{}_AAAAAA_L00{}_R{}_001.fastq.gz'.format(sample_name,info[-3],'2')
                read_3_flag = 1
            elif info[-2] == '2' and read_3_flag == 0:
                new_file_name = '{}_AAAAAA_L00{}_R{}_001.fastq.gz'.format(sample_name,info[-3],'2')
            else:
                raise Exception, 'File {} could not be processed'.format(sequence_file)

            seqtk_gzip_cmd = ('/nfs/goldstein/software/seqtk/seqtk seq -Q64 -V {}'
                                ' | /nfs/goldstein/software/pigz-2.3.1/pigz -p 4 > {}/{}').format(sequence_file,dest,new_file_name)
            print seqtk_gzip_cmd
            if run_script:
                subprocess.check_call(seqtk_gzip_cmd,shell=True)

            rm_cmd = 'rm {}'.format(sequence_file)
            runCmd(rm_cmd)

        counter += 1

def runCmd(cmd):
    print cmd
    if run_script:
        subprocess.check_call(shlex.split(cmd))

if __name__ == '__main__':
    main()
