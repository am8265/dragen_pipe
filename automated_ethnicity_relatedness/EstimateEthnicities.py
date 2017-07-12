#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Automate calculating ethnicity probabilities for samples with production
pipeline
"""

import argparse
import re
import os
import sys
sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "import", "src"))
import ProcessSamples
import logging
from logging import handlers
from ethnicity_db_statements import *
from waldb_globals import *

cfg = get_cfg()

def find_vcf(data_directory, sample_name):
    """Try to find a sample's VCF (gzipped or not), otherwise return None
    """
    prefix = os.path.join(data_directory, "combined", sample_name)
    for vcf_name in (".analysisReady.vcf", ".analysisReady.annotated.vcf"):
        for suffix in ("", ".gz"):
            vcf = prefix + vcf_name + suffix
            if os.path.isfile(vcf):
                return vcf
    return None

class CalculateEthnicityProbabilities(ProcessSamples.ProcessSamples):
    def __init__(self, output_directory, max_samples_concurrently=200, **kwargs):
        super(CalculateEthnicityProbabilities, self).__init__(
            max_samples_concurrently=max_samples_concurrently, **kwargs)
        self.output_directory = output_directory
        if not os.path.isdir(output_directory):
            os.makedirs(output_directory)
        self.pattern = (
            r"^(?:CreatePed|CalculateProbabilities)\.{sample_name}\.{0}.{1}$")

    def _get_samples(self):
        db = get_connection("seqdb")
        try:
            samples = []
            cur = db.cursor()
            query = PROD_GET_SAMPLES
            cur.execute(query)
            for row in cur.fetchall():
                samples.append((row[0], row[1], row[2:]))
            return samples
        finally:
            if db.open:
                db.close()

    def _get_command(self, sample_name, *args):
        prep_id, sample_type, data_directory = args
        vcf = find_vcf(data_directory, sample_name)
        bam = os.path.join(
            data_directory, "combined", "combined_rmdup_realn_recal.bam")
        if vcf and os.path.isfile(bam):
            cmd = ("luigi --module estimate_ethnicities CalculateProbabilities "
                   "--sample-name {sample_name} --prep-id {prep_id} --sample-type "
                   "{sample_type} --vcf {vcf} --bam {bam} --output-directory "
                   "{output_directory}/{{sample_name}}.{{prep_id}}.{{sample_type}} "
                   "--logging-conf-file logging.conf{run_locally}".format(
                       sample_name=sample_name, prep_id=prep_id,
                       sample_type=sample_type, vcf=vcf, bam=bam,
                       output_directory=self.output_directory,
                       run_locally=" --run-locally" if self.run_locally else ""))
        else:
            logger = logging.getLogger("ProcessSamples")
            logger.debug("Couldn't find VCF/BAM for {sample_name}:{sample_type}:{prep_id}".
                         format(sample_name=sample_name,
                                sample_type=sample_type, prep_id=prep_id))
            cmd = "/bin/false"
        timeout = 7200 if sample_type == "Exome" else 14400
        return cmd, timeout

def main(num_tasks, run_locally, output_directory):
    formatter = logging.Formatter(cfg.get("logging", "format"))
    logger = logging.getLogger("ProcessSamples")
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.INFO)
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    C = CalculateEthnicityProbabilities(
        output_directory, max_samples_concurrently=num_tasks, run_locally=run_locally)
    C.process_samples()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=CustomFormatter)
    parser.add_argument("-n", "--num_tasks", default=200, type=int,
                        help="the number of samples to process simultaneously")
    parser.add_argument("-l", "--run_locally", default=False,
                        action="store_true",
                        help="run locally instead of on the cluster")
    parser.add_argument("-o", "--output_directory",
                        default="/nfs/seqscratch_ssd/bc2675/annodb_ethnicities",
                        help="the directory to output data to")
    args = parser.parse_args()
    main(args.num_tasks, args.run_locally, args.output_directory)
