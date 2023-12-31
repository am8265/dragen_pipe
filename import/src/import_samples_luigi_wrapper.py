#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Automate import of samples to the DRAGEN DB
"""

import argparse
import re
import os
import ProcessSamples
import logging
import sys
from logging import handlers
from waldb_globals import *
from db_statements import *

cfg = get_cfg()

class ImportSamples(ProcessSamples.ProcessSamples):
    def __init__(self, database="waldb6", seqscratch="_ssd",
                 force_failed_samples=False,
                 timeout_exome=1800, timeout_genome=7200,
                 timeout_custom_capture=1800, sample_names=None,
                 run_locally=False, workers=40, qdel_jobs=True,
                 local_scheduler=False, **kwargs):
        super(ImportSamples, self).__init__(
            max_samples_concurrently=1, ignore_priority=False,
            force_failed_samples=force_failed_samples,
            sample_names=sample_names,
            run_locally=run_locally, qdel_jobs=qdel_jobs,
            local_scheduler=local_scheduler, **kwargs)
        self.pattern = (
            r"^(?:ParseVCF|LoadGQData|LoadDPData)\.{sample_name}\."
            r"{0}\.chr(?:1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|"
            r"16|17|18|19|20|21|22|X|Y|MT)$")
        self.database = database
        self.seqscratch = seqscratch
        self.workers = workers
        self.timeout = {"exome":timeout_exome, "genome":timeout_genome,
                        "custom_capture":timeout_custom_capture}

    def _get_samples(self):
        db = get_connection(self.database)
        try:
            cur = db.cursor()
            query = GET_SAMPLES_TO_IMPORT.format(
                failed_samples_clause="" if self.force_failed_samples
                else " AND sample_failure = 0",
            sample_name_clause=" AND sample_name LIKE "
                "'{sample_names}'".format(sample_names=self.sample_names)
                if self.sample_names else "")
            cur.execute(query)
            samples = []
            for row in cur.fetchall():
                samples.append((row[0], row[1], row[2], row[3:]))
        finally:
            if db.open:
                db.close()
        return samples

    def _get_command(self, pseudo_prepid, sample_name, *args):
        sample_id, sample_type = args
        cmd = ("luigi --module data_import_pipeline ImportSample --sample-id "
               "{sample_id} --seqscratch {seqscratch} --workers {workers} "
               "--database {database} "
               "--dont-remove-tmp-dir-if-failure "
               "--logging-conf-file {logging_conf_file}{run_locally}"
               "{local_scheduler}".format(
                   sample_id=sample_id, seqscratch=self.seqscratch,
                   workers=self.workers,
                   logging_conf_file="logging.conf", database=self.database,
                   run_locally=" --run-locally" if self.run_locally else "",
                   local_scheduler=" --local-scheduler" if self.local_scheduler
                   else ""))
        #return cmd, self.timeout[sample_type]

        ######## really, why the heck would you need to be 'careful' just block it if it's dangerous?!?
        return cmd, None # we have to be very careful to not ever load multiple
                         # samples, or multiple variants can be assigned the same
                         # variant_id

def main(database, seqscratch, force_failed_samples, sample_names,
         run_locally, workers, local_scheduler, debug_level,
         timeout_exome, timeout_genome, timeout_custom_capture):
    print("i'm deprecated")
    exit(1)
    # set up logging
    formatter = logging.Formatter(cfg.get("logging", "format"))
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    root_handler = handlers.TimedRotatingFileHandler(
        "/nfs/seqscratch_ssd/informatics/dragen_import/dragen_import_log.txt",
        when="midnight", backupCount=10)
    root_handler.setLevel(logging.DEBUG)
    root_handler.setFormatter(formatter)
    root_logger.addHandler(root_handler)
    logger = logging.getLogger("ProcessSamples")
    logger.setLevel(debug_level)
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(debug_level)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    # import samples
    import_samples = ImportSamples(
        database=database, seqscratch=seqscratch,
        force_failed_samples=force_failed_samples,
        timeout_exome=timeout_exome, timeout_genome=timeout_genome,
        timeout_custom_capture=timeout_custom_capture,
        sample_names=sample_names,
        qdel_jobs=not run_locally, run_locally=run_locally,
        workers=workers, local_scheduler=local_scheduler)
    import_samples.process_samples()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=CustomFormatter)
    parser.add_argument("-d", "--database", default="waldb6",
                        choices=["waldb", "waldb2", "waldb4", "waldb1", "waldb6"],
                        help="the database to import to")
    parser.add_argument("-s", "--seqscratch",
                        choices=["09", "10", "11", "_ssd"], default="_ssd",
                        help="the scratch drive to use")
    parser.add_argument("--force_failed_samples", default=False,
                        action="store_true",
                        help="try to import a sample even if it failed in a "
                        "previous run (otherwise such samples are ignored)")
    parser.add_argument("--sample_names", help="only attempt to import samples "
                        "with a sample_name matching the format string")
    parser.add_argument("--run_locally", default=False, action="store_true",
                        help="run locally instead of on the cluster")
    parser.add_argument("-w", "--workers", type=int, default=40,
                        help="the number of luigi workers to use (limit to "
                        "restrict DB connections)")
    parser.add_argument("--local_scheduler", default=False, action="store_true",
                        help="use the local luigi scheduler instead of "
                        "luigi daemon")
    parser.add_argument("--timeout_exome", default=1800, type=int,
                        help="the timeout value for exomes")
    parser.add_argument("--timeout_genome", default=7200, type=int,
                        help="the timeout value for genomes")
    parser.add_argument("--timeout_custom_capture", default=1800, type=int,
                        help="the timeout value for custom capture samples")
    parser.add_argument("--level", default="INFO", choices=LOGGING_LEVELS,
                        action=DereferenceKeyAction,
                        help="specify the logging level to use")
    args = parser.parse_args()
    main(args.database, args.seqscratch, args.force_failed_samples,
         args.sample_names, args.run_locally, args.workers,
         args.local_scheduler, args.level, args.timeout_exome,
         args.timeout_genome, args.timeout_custom_capture)
