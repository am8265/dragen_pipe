#!/usr/bin/env python
"""
Delete sample(s)' variant calls from WalDB and reset sample_finished flag
"""

import multiprocessing
import argparse
import logging
import logging
import sys
from waldb_globals import *

cfg = get_cfg()
CHROMs = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
          "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"]

def get_sample_ids(input_sample_ids, sample_name_string):
    """return all sample_id(s) found
    """
    sample_ids = []
    db = get_connection("waldb4")
    try:
        cur = db.cursor()
        if input_sample_ids:
            for sample_id in input_sample_ids:
                cur.execute("SELECT sample_id FROM sample WHERE sample_id = {}".
                            format(sample_id))
                if cur.fetchone():
                    sample_ids.append(sample_id)
        if sample_name_string:
            cur.execute("SELECT sample_id FROM sample WHERE sample_name LIKE "
                        "'{}'".format(sample_name_string))
            for row in cur.fetchall():
                sample_ids.append(row[0])
        return sample_ids
    finally:
        if db.open:
            db.close()

@timer(sys.stderr)
def delete_sample_chromosome(sample_id, CHROM):
    """delete a given sample_id's calls on the specified chromosome
    """
    db = get_connection("waldb4")
    seqdb = get_connection("seqdb")
    try:
        cur = db.cursor()
        seq_cur = seqdb.cursor()
        cur.execute("SELECT prep_id FROM sample WHERE sample_id = {}".format(
            sample_id))
        prep_id = cur.fetchone()[0]
        seq_cur.execute("SELECT id FROM dragen_pipeline_step_desc WHERE "
                        "step_name = 'Imported Chromosome {} Variant Data'".
                        format(CHROM))
        pipeline_step_id = seq_cur.fetchone()[0]
        cur.execute("DELETE FROM called_variant_chr{} WHERE sample_id = {}".
                    format(CHROM, sample_id))
        cur.execute("UPDATE sample SET sample_finished = 0 WHERE sample_id = {}".
                    format(sample_id))
        seq_cur.execute("DELETE FROM dragen_pipeline_step WHERE pseudo_prepid = {} "
                        "AND pipeline_step_id = {}".format(
                            prep_id, pipeline_step_id))
        db.commit()
        seqdb.commit()
    finally:
        if db.open:
            db.close()
        if seqdb.open:
            seqdb.close()

@timer(sys.stderr)
def delete_samples_chromosome(sample_ids, CHROM):
    """delete all samples' calls on the specified chromosome
    """
    db = get_connection("waldb4")
    seqdb = get_connection("seqdb")
    try:
        cur = db.cursor()
        seq_cur = seqdb.cursor()
        prep_ids = []
        for sample_id in sample_ids:
            cur.execute("SELECT prep_id FROM sample WHERE sample_id = {}".format(
                sample_id))
            prep_ids.append(str(cur.fetchone()[0]))
        seq_cur.execute("SELECT id FROM dragen_pipeline_step_desc WHERE "
                        "step_name = 'Imported Chromosome {} Variant Data'".
                        format(CHROM))
        pipeline_step_id = seq_cur.fetchone()[0]
        sample_ids_str = ",".join([str(sample_id) for sample_id in sample_ids])
        cur.execute("DELETE FROM called_variant_chr{} WHERE sample_id IN ({})".
                    format(CHROM, sample_ids_str))
        cur.execute("UPDATE sample SET sample_finished = 0 WHERE sample_id IN "
                    "({})".format(sample_ids_str))
        seq_cur.execute("DELETE FROM dragen_pipeline_step WHERE "
                        "pipeline_step_id IN (0, {}) AND pseudo_prepid IN ({})".format(
                            pipeline_step_id, ",".join(prep_ids)))
        db.commit()
        seqdb.commit()
    finally:
        if db.open:
            db.close()
        if seqdb.open:
            seqdb.close()

@timer(sys.stderr)
def delete_sample(sample_id):
    """delete a given sample's calls, reset sample_finished = 0, and delete the
    corresponding record in SequenceDB
    """
    pool = multiprocessing.Pool()
    for CHROM in CHROMs:
        pool.apply_async(delete_sample_chromosome, args=(sample_id, CHROM))
    pool.close()
    pool.join()

@timer(sys.stderr)
def delete_samples(sample_ids):
    """delete all samples' calls, reset sample_finished = 0, and delete the
    corresponding record in SequenceDB
    """
    pool = multiprocessing.Pool()
    for CHROM in CHROMs:
        pool.apply_async(delete_samples_chromosome, args=(sample_ids, CHROM))
    pool.close()
    pool.join()

def main(sample_ids):
    """delete all specified sample_ids' calls, delete the pipeline step entry in
    SeqDB, and reset sample_finished = 0 in sample table
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(cfg.get("logging", "format"))
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    nsamples = len(sample_ids)
    delete_samples(sample_ids)
    #for x, sample_id in enumerate(sample_ids):
    #    logger.info("starting sample #{x}/{n} ({sample_id})".format(
    #        x=x + 1, n=nsamples, sample_id=sample_id))
    #    delete_sample(sample_id)
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=CustomFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sample_id", type=int, action="append",
                       help="a single sample_id to delete")
    group.add_argument("--sample_name_string", help="a string to match the "
                       "sample_name field against, e.g. 'SRR%'")
    args = parser.parse_args()
    main(get_sample_ids(args.sample_id, args.sample_name_string))
