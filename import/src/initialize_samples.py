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
from db_statements import (
    GET_SAMPLE_ID, GET_SAMPLE_PREPID, INSERT_SAMPLE, GET_PREPID)

def initialize_samples(samples_fh, logging_level):
    """Take the list of samples in the fh and initialize all samples in the
    dragen db
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging_level)
    seqdb = get_connection("seqdb")
    dragen = get_connection("dragen")
    try:
        seqdb_cur = seqdb.cursor()
        dragen_cur = dragen.cursor()
        samples_metadata = []
        for line in samples_fh:
            (family_id, sample_name, paternal_id, maternal_id, sex, phenotype,
             sample_type, capture_kit) = line.strip().split("\t")
            sample_metadata = {
                "sample_name":sample_name, "sample_type":sample_type.lower(),
                "capture_kit":capture_kit}
            dragen_cur.execute(GET_SAMPLE_ID.format(**sample_metadata))
            if dragen_cur.fetchone():
                logging.debug("{sample_name} is in DB".format(**sample_metadata))
                continue
            query = GET_SAMPLE_PREPID.format(
                sample_name=sample_name, sample_type=sample_type,
                capture_kit=capture_kit)
            seqdb_cur.execute(query)
            rows = seqdb_cur.fetchall()
            if len(rows) == 1:
                pseudo_prep_id, prep_id = rows[0]
                if pseudo_prep_id is None:
                    query_two = GET_PREPID.format(prepid=prep_id)
                    seqdb_cur.execute(query_two)
                    row = seqdb_cur.fetchone()
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
                samples_metadata["prep_id"] = pseudo_prep_id
                logging.info("inserting {sample_name}".format(**sample_metadata))
                dragen_cur.execute(INSERT_SAMPLE.format(**sample_metadata))
                dragen.commit()
            else:
                logging.warning("Found {} records with query:{}".format(
                    len(rows), query))
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
