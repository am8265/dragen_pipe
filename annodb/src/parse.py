#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Take an annotated VCF, insert novel variants, and create output tables to
subsequently load into the database

Brett Copeland <bc2675@cumc.columbia.edu>
3/2016
"""

import argparse
import sys
import MySQLdb
import subprocess
import shlex
from ConfigParser import SafeConfigParser
from utilities import *

CHROMs = set([str(chrom) for chrom in range(1, 23)] + ["X", "Y", "MT"])
SAMPLE_EXISTS_QUERY = """
SELECT sample_name
FROM sample
WHERE sample_id = {sample_id}
"""
CNF = "/nfs/goldstein/software/dragen/dragen.cnf"
CHROM_FH_FORMAT = "{output_directory}/{sample_name}_chr{CHROM}.vcf"
config_parser = SafeConfigParser()
config_parser.read(CNF)
script_path = get_script_path()
PARSE_CHROM_CMD = (
    "{script_path}/parse_vcf_chromosome.py {VCF} "
    "{output_directory} {sample_id} {CHROM}")

def submit_chromosome(sample_name, sample_id, output_directory, CHROM):
    """submit a job to the cluster to process the specified chromosome
    """
    full_cmd = ("qsub -S /bin/sh -N {sample_name}.parse_vcf_chr{CHROM} -o "
                "{output_directory}/{sample_name}.parse_vcf_chr{CHROM}.out -e "
                "{output_directory}/{sample_name}.parse_vcf_chr{CHROM}.err -wd "
                "{output_directory} -b y " + PARSE_CHROM_CMD).format(
                    CHROM=CHROM, sample_name=sample_name, sample_id=sample_id,
                    output_directory=output_directory,
                    script_path=script_path, VCF=CHROM_FH_FORMAT.format(
                        output_directory=output_directory,
                        sample_name=sample_name, CHROM=CHROM))
    print(full_cmd)
    p = subprocess.Popen(shlex.split(full_cmd), stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    rc = p.returncode
    print("stdout={}".format(stdout))
    print("stderr={}".format(stderr))
    print("rc={}".format(rc))

def main(input_vcf_fh, output_directory, sample_id, log_fh):
    """split the VCF by chromosome and submit jobs to the cluster to parse each
    chromosome independently (concurrently)
    """
    db = MySQLdb.connect(
        read_default_file=config_parser.get("annodb", "DB_CONFIG_FILE"),
        read_default_group="clientdragen01")
    cur = db.cursor()
    cur.execute(SAMPLE_EXISTS_QUERY.format(sample_id=sample_id))
    row = cur.fetchone()
    db.close()
    if row:
        sample_name = row[0]
    else:
        log_fh.write("error parsing VCF: sample_id {sample_id} "
                     "does not exist\n".format(sample_id=sample_id))
    current_CHROM = -1
    for line in input_vcf_fh:
        if line.startswith("#CHROM"):
            break
    for line in input_vcf_fh:
        CHROM = line.split("\t")[0]
        if CHROM != current_CHROM:
            if current_CHROM in CHROMs:
                chrom_fh.close()
                submit_chromosome(
                    sample_name, sample_id, output_directory, current_CHROM)
            current_CHROM = CHROM
            if current_CHROM in CHROMs:
                chrom_fh = open(CHROM_FH_FORMAT.format(
                    output_directory=output_directory, sample_name=sample_name,
                    CHROM=CHROM), "w")
            else:
                break
        chrom_fh.write(line)
    input_vcf_fh.close()
    if not chrom_fh.closed:
        chrom_fh.close()
    if log_fh is not sys.stderr:
        log_fh.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=CustomFormatter, description=__doc__)
    parser.add_argument("INPUT_VCF", type=fh, help="the VCF to parse")
    parser.add_argument("OUTPUT_DIRECTORY", type=make_directory_if_not_exists,
                        help="the directory to output to")
    parser.add_argument("SAMPLE_ID", type=natural_number,
                        help="the sample_id of the sample whose VCF this is")
    parser.add_argument("--log", type=argparse.FileType("a"),
                        default=sys.stderr,
                        help="the log file to write errors to")
    args = parser.parse_args()
    main(args.INPUT_VCF, args.OUTPUT_DIRECTORY, args.SAMPLE_ID, args.log)
