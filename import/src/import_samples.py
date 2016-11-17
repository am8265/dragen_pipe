#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Automate import of samples to the DRAGEN DB
"""

import argparse
import re
import os
import ProcessSamples
from dragen_globals import *
from db_statements import *

class ImportSamples(ProcessSamples.ProcessSamples):
    def __init__(self, qdel_jobs=True, **kwargs):
        super(ImportSamples, self).__init__(
            max_samples_concurrently=1, qdel_jobs=qdel_jobs, **kwargs)
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
        cmd = ("./data_import_pipeline.py ImportSample --sample-id {sample_id} "
               "--workers 75 --logging-conf-file {logging_conf_file}{run_locally}".
               format(
                   sample_id=sample_id, logging_conf_file="logging.conf",
                   run_locally=" --run-locally" if self.kwargs["run_locally"] else ""))
        if sample_type in ("exome", "custom_capture"):
            # allow 10 minutes for exomes
            timeout = 600
        elif sample_type in ("genome", "merged"):
            # 2 hours for genomes
            timeout = 7200
        return cmd, timeout

def main(run_locally):
    import_samples = ImportSamples(
        qdel_jobs=not run_locally, run_locally=run_locally,
        stdout=os.devnull)
    import_samples.process_samples()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=CustomFormatter)
    parser.add_argument("--run_locally", default=False, action="store_true",
                        help="run locally instead of on the cluster")
    args = parser.parse_args()
    main(args.run_locally)
