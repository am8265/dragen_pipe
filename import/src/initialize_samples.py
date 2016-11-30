#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Initialize a list of samples in the DRAGEN db
"""

import argparse
import MySQLdb
import sys
import os
import logging
from dragen_globals import *
from db_statements import *

def initialize_samples(samples_fh, logging_level):
    """Take the list of samples in the fh and initialize all samples in the
    dragen db
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging_level)
    seqdb = get_connection("seqdb")
    dragen = get_connection("dragen")
    try:
        seq_cur = seqdb.cursor()
        query = GET_PIPELINE_SAMPLE_INITIALIZED_ID
        seq_cur.execute(query)
        pipeline_step_id = seq_cur.fetchone()[0]
        dragen_cur = dragen.cursor()
        for line in samples_fh:
            (family_id, sample_name, paternal_id, maternal_id, sex, phenotype,
             sample_type, capture_kit) = line.strip().split("\t")
            sample_metadata = {
                "sample_name":sample_name, "sample_type":sample_type.lower(),
                "capture_kit":capture_kit}
            query = GET_SAMPLE_ID.format(**sample_metadata)
            dragen_cur.execute(query)
            if dragen_cur.fetchone():
                logging.debug("{sample_name} is in DB".format(**sample_metadata))
                continue
            query = GET_SAMPLE_PREPID.format(
                sample_name=sample_name, sample_type=sample_type,
                capture_kit=capture_kit)
            seq_cur.execute(query)
            rows = seq_cur.fetchall()
            if len(rows) == 1:
                pseudo_prep_id, prep_id, priority = rows[0]
                if pseudo_prep_id is None:
                    query = GET_PREPID.format(prepid=prep_id)
                    seq_cur.execute(query)
                    row = seq_cur.fetchone()
                    if row:
                        pseudo_prep_id = row[0]
                data_directory = (
                    "/nfs/fastq16/ALIGNMENT/BUILD37/DRAGEN/{sample_type}/"
                    "{sample_name}.{pseudo_prep_id}".format(
                        sample_type=sample_type.upper(),
                        sample_name=sample_name, pseudo_prep_id=pseudo_prep_id))
                if not os.path.isdir(data_directory):
                    logging.warning("couldn't find data directory {}".format(
                        data_directory))
                    continue
                sample_metadata["prep_id"] = pseudo_prep_id
                sample_metadata["priority"] = priority
                logging.info("inserting {sample_name}".format(**sample_metadata))
                query = INITIALIZE_SAMPLE_PIPELINE_STEP.format(
                    pipeline_step_id=pipeline_step_id, **sample_metadata)
                seq_cur.execute(query)
                query = INSERT_SAMPLE.format(**sample_metadata)
                dragen_cur.execute(query)
                seqdb.commit()
                dragen.commit()
            else:
                logging.warning("Found {} records with query:{}".format(
                    len(rows), query))
    except:
        logging.error("Query failed:\n" + query)
        raise
    finally:
        samples_fh.close()
        if seqdb.open:
            seqdb.close()
        if dragen.open:
            dragen.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=CustomFormatter)
    parser.add_argument("SAMPLES_FN", type=argparse.FileType("r"),
                        help="the sample list file")
    parser.add_argument("--level", default="DEBUG", action=DereferenceKeyAction,
                        choices=LOGGING_LEVELS, help="the logging level to use")
    args = parser.parse_args()
    initialize_samples(args.SAMPLES_FN, args.level)
