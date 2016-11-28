#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Automate import of samples to the DRAGEN DB
"""

import argparse
import re
import os
import ProcessSamples
import logging
from logging import handlers
from dragen_globals import *
from db_statements import *

cfg = get_cfg()
formatter = logging.Formatter(cfg.get("logging", "format"))
root_logger = logging.getLogger()
root_logger.setLevel(logging.DEBUG)
root_handler = handlers.TimedRotatingFileHandler(
    "/nfs/seqscratch09/dragen_import/dragen_import_log.txt", when="midnight",
    backupCount=10)
root_handler.setLevel(logging.DEBUG)
root_handler.setFormatter(formatter)
root_logger.addHandler(root_handler)
logger = logging.getLogger("ProcessSamples")
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stderr)
handler.setLevel(logging.DEBUG)
handler.setFormatter(formatter)
logger.addHandler(handler)

class ImportSamples(ProcessSamples.ProcessSamples):
    def __init__(self, run_locally=False, qdel_jobs=True,
                 local_scheduler=False, **kwargs):
        super(ImportSamples, self).__init__(
            max_samples_concurrently=1, run_locally=run_locally,
            qdel_jobs=qdel_jobs, local_scheduler=local_scheduler, **kwargs)
        self.pattern = (
            r"^(?:ParseVCF|LoadGQData|LoadDPData)\.{sample_name}\."
            r"{0}\.chr(?:1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|"
            r"16|17|18|19|20|21|22|X|Y|MT)$")

    def _get_samples(self):
        db = get_connection("dragen")
        try:
            cur = db.cursor()
            cur.execute(GET_SAMPLES_TO_IMPORT)
            samples = []
            for row in cur.fetchall():
                samples.append((row[0], row[1], row[2:]))
        finally:
            if db.open:
                db.close()
        return samples

    def _get_command(self, sample_name, *args):
        sample_id, sample_type = args
        cmd = ("luigi --module data_import_pipeline ImportSample --sample-id "
               "{sample_id} --workers 75 --dont-remove-tmp-dir-if-failure "
               "--logging-conf-file {logging_conf_file}{run_locally}"
               "{local_scheduler}".format(
                   sample_id=sample_id, logging_conf_file="logging.conf",
                   run_locally=" --run-locally" if self.run_locally else "",
                   local_scheduler=" --local-scheduler" if self.local_scheduler
                   else ""))
        if sample_type in ("exome", "custom_capture"):
            # allow 10 minutes for exomes
            timeout = 600
        elif sample_type in ("genome", "merged"):
            # 2 hours for genomes
            timeout = 7200
        return cmd, timeout

def main(run_locally, local_scheduler):
    import_samples = ImportSamples(
        qdel_jobs=not run_locally, run_locally=run_locally,
        local_scheduler=local_scheduler)
    import_samples.process_samples()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=CustomFormatter)
    parser.add_argument("--run_locally", default=False, action="store_true",
                        help="run locally instead of on the cluster")
    parser.add_argument("--local_scheduler", default=False, action="store_true",
                        help="use the local luigi scheduler instead of "
                        "luigi daemon")
    args = parser.parse_args()
    main(args.run_locally, args.local_scheduler)
