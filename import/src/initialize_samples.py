#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Initialize a list of samples in the DRAGEN db
"""

import argparse
import MySQLdb
import sys
import os
from dragen_globals import *
from db_statements import (
    GET_SAMPLE_ID, GET_SAMPLE_PREPID, INSERT_SAMPLE, GET_PREPID)

def initialize_samples(samples_fh):
    """Take the list of samples in the fh and initialize all samples in the
    dragen db
    """
    seqdb = get_connection("seqdb")
    dragen = get_connection("dragen")
    try:
        seqdb_cur = seqdb.cursor()
        dragen_cur = dragen.cursor()
        errors = []
        samples_metadata = []
        for line in samples_fh:
            (family_id, sample_name, paternal_id, maternal_id, sex, phenotype,
             sample_type, capture_kit) = line.strip().split("\t")
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
                if not os.path.isdir(
                    "/nfs/fastq16/ALIGNMENT/BUILD37/DRAGEN/{sample_type}/"
                    "{sample_name}.{pseudo_prep_id}".format(
                        sample_type=sample_type.upper(),
                        sample_name=sample_name, pseudo_prep_id=pseudo_prep_id)):
                    errors.append((query, -1))
                    continue
                samples_metadata.append({
                    "sample_name":sample_name,
                    "sample_type":sample_type.lower(), "capture_kit":capture_kit,
                    "prep_id":pseudo_prep_id})
            else:
                errors.append((query, len(rows)))
        if errors:
            for query, nrecords in errors:
                print("error with query ({} records):{}".format(
                    nrecords, query))
            sys.exit(1)
        for sample_metadata in samples_metadata:
            dragen_cur.execute(GET_SAMPLE_ID.format(**sample_metadata))
            if not dragen_cur.fetchone():
                dragen_cur.execute(INSERT_SAMPLE.format(**sample_metadata))
        dragen.commit()
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
    args = parser.parse_args()
    initialize_samples(args.SAMPLES_FN)
