#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Automate running the GATK pipeline, annotation, binning, QC, etc.
defined in the gatk_pipe luigi script
"""

import argparse
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.realpath(__file__)))), "import", "src"))
import logging
from time import time as curr_time
from logging import handlers
from ProcessSamples import ProcessSamples
from dragen_db_statements import *
from waldb_globals import *

cfg = get_cfg()

class RunGATK(ProcessSamples):
    def __init__(self, max_samples_concurrently, workers,
                 additional_sample_requirements, no_retry, luigi_args):
        super(RunGATK, self).__init__(max_samples_concurrently)
        self.workers = workers
        query = GET_SAMPLES
        self.query = (query + additional_sample_requirements if
                      additional_sample_requirements else query)
        self.luigi_args = " " + " ".join(luigi_args) if luigi_args else ""
        self.no_retry = no_retry
        self.last_query_time = None
        self.time_between_queries = 60
        seqdb = get_connection("seqdb")
        try:
            seq_cur = seqdb.cursor()
            seq_cur.execute(GET_STEP_NAMES)
            self.pattern = (r"^(?:" +
                            "|".join([row[0] for row in seq_cur.fetchall()]) +
                            r")\.{sample_name}\.{pseudo_prepid}$")
        finally:
            if seqdb.open:
                seqdb.close()
        self.qdel_jobs = True

    def _get_samples(self):
        # execute the GET_SAMPLES query, but not more than every X seconds
        current_time = curr_time()
        samples = []
        if (self.last_query_time and (current_time - self.last_query_time) <
            self.time_between_queries):
            return samples
        self.last_query_time = current_time
        seqdb = get_connection("seqdb")
        try:
            seq_cur = seqdb.cursor()
            seq_cur.execute(self.query)
            rows = seq_cur.fetchall()
            for row in rows:
                if self.no_retry:
                    pseudo_prepid = row[0]
                    seq_cur.execute(ANY_STEP_FAILED.format(
                        pseudo_prepid=pseudo_prepid))
                    # only run the sample if none of the steps failed
                    run_sample = not bool(seq_cur.fetchone())
                else:
                    run_sample = True
                if run_sample:
                    samples.append((row[0], row[1], row[2], ()))
            return samples
        finally:
            if seqdb.open:
                seqdb.close()

    def _get_command(self, pseudo_prepid, sample_name):
        """pass on a minimal set of parameters, but any additional arbitrary ones
        specified in luigi_args
        """
        return ("luigi --module gatk_pipe ArchiveSample --logging-conf-file "
                "logging.conf --dont-remove-tmp-dir-if-failure --pseudo-prepid "
                "{pseudo_prepid} --workers {workers}{luigi_args}".format(
                    pseudo_prepid=pseudo_prepid, workers=self.workers,
                    luigi_args=self.luigi_args), None)

def run_gatk(max_samples_concurrently, workers, additional_sample_requirements,
             no_retry, debug_level, luigi_args):
    """Run the gatk_pipe.py code with the parameters specified
    @max_samples - the max number of samples to process simultaneously
    @workers - the number of workers each gatk_pipe run should have
    @additional_sample_requirements - any additions to the SQL query
    @luigi_args - any additional parameters to be passed on to luigi
    """
    # Allow group write permissions for all created files
    os.umask(0002)
    formatter = logging.Formatter(cfg.get("logging", "format"))
    root_logger = logging.getLogger()
    root_logger.setLevel(debug_level)
    root_handler = handlers.TimedRotatingFileHandler(
        "/nfs/seqscratch09/dragen_gatk/gatk_pipe.txt",
        when="midnight", backupCount=10)
    root_handler.setLevel(debug_level)
    root_handler.setFormatter(formatter)
    root_logger.addHandler(root_handler)
    logger = logging.getLogger("ProcessSamples")
    logger.setLevel(debug_level)
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(debug_level)
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    gatk = RunGATK(
        max_samples_concurrently, workers,
        additional_sample_requirements, no_retry, luigi_args)
    gatk.process_samples()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=CustomFormatter, description=__doc__)
    parser.add_argument("-m", "--max_samples_concurrently", type=int, default=50,
                        help="process at most this many samples at once")
    parser.add_argument("-w", "--workers", type=int, default=2,
                        help="use this many workers per sample")
    parser.add_argument("--additional_sample_requirements", help="add additional "
                        "SQL query logic to restrict samples processed")
    parser.add_argument("-n", "--no_retry", default=False, action="store_true",
                        help="don't attempt to retry a sample if any of its "
                        "steps have failed previously")
    parser.add_argument("--level", default="INFO", choices=LOGGING_LEVELS,
                        action=DereferenceKeyAction,
                        help="specify the logging level to use")
    args, luigi_args = parser.parse_known_args()
    run_gatk(args.max_samples_concurrently, args.workers,
             args.additional_sample_requirements,
             args.no_retry, args.level, luigi_args)
