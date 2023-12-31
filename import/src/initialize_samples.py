#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Initialize samples that have finished the DRAGEN pipeline in WalDB
"""

# just run INITIALIZE_SAMPLE & INITIALIZE_SAMPLE_SEQDB over GET_SAMPLES_TO_INITIALIZE - just a self join of dps for samples with 31 and not 32...?!?
# then 'replaces into' dps as initialised... and inserts into sample...?!?

import argparse
import MySQLdb
import logging
from sys import stderr
from waldb_globals import *
from db_statements import *

cfg = get_cfg()

@timer(stderr)
def initialize_samples(database, level=logging.DEBUG):
    logger = logging.getLogger(__name__)
    logger.setLevel(level)
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(level)
    formatter = logging.Formatter(cfg.get("logging", "format"))
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    seqdb = get_connection("seqdb")
    db = get_connection(database)
    initialized_count = 0
    try:

        seq_cur = seqdb.cursor()
        cur = db.cursor()

        #### why separate queries
        seq_cur.execute(                GET_PIPELINE_STEP_ID.format(    step_name="ArchiveSample")      )
        sample_archived_step_id     =   seq_cur.fetchone()[0]
        seq_cur.execute(                GET_PIPELINE_STEP_ID.format(    step_name="Sample Initialized in DB")   )
        sample_initialized_step_id  =   seq_cur.fetchone()[0]

        seq_cur.execute(GET_SAMPLES_TO_INITIALIZE.format( sample_archived_step_id=sample_archived_step_id,
            sample_initialized_step_id=sample_initialized_step_id)
        )
        prep_ids = [row[0] for row in seq_cur.fetchall()]
        nsamples = len(prep_ids)

        logger.info("Found {nsamples} samples to initialize".format(nsamples=nsamples))
        for prep_id in prep_ids:

            if (111000 <= prep_id < 120000) or (prep_id >= 838600):
                logger.debug("Ignoring test sample: {}".format(prep_id))
                continue

            seq_cur.execute(GET_SAMPLE_METADATA.format(prep_id=prep_id))
            rows = seq_cur.fetchall()

            if not rows:
                logger.warning("Couldn't get metadata for {prep_id}".format(
                    prep_id=prep_id))
            elif len(rows) > 1:
                logger.warning("Found multiple records for {prep_id}".format(
                    prep_id=prep_id))
            else:
                sample_name, sample_type, capture_kit, scratch, priority = rows[0]

                if sample_name.startswith("pgm"):
                    logger.info( "Ignoring {sample_name}:{sample_type}:" "{capture_kit}:{prep_id}".format(
                      sample_name=sample_name, sample_type=sample_type, capture_kit=capture_kit, prep_id=prep_id) )
                    continue

                if sample_type.lower() == "genome":
                    capture_kit = "N/A"
                logger.debug( "Initializing {sample_name}:{sample_type}:{capture_kit}:{prep_id}".format(
                  sample_name=sample_name, sample_type=sample_type, capture_kit=capture_kit, prep_id=prep_id) )

                try:
                    cur.execute(INITIALIZE_SAMPLE.format( sample_name=sample_name, sample_type=sample_type,
                      capture_kit=capture_kit, prep_id=prep_id, priority=priority))
                    add_to_count = True
                except MySQLdb.IntegrityError:
                    add_to_count = False
                    logger.warning("Prep ID {prep_id} is already initialized".format( prep_id=prep_id) )
                version = get_pipeline_version()
                seq_cur.execute(INITIALIZE_SAMPLE_SEQDB.format( pseudo_prepid=prep_id, version=version,
                  sample_initialized_step_id=sample_initialized_step_id))
                if add_to_count:
                    initialized_count += 1
        db.commit()
        seqdb.commit()
    except:
        db.rollback()
        seqdb.rollback()
        raise
    finally:
        if seqdb.open:
            seqdb.close()
        if db.open:
            db.close()
        logger.info("Initialized {initialized_count}/{nsamples} samples".format(
            initialized_count=initialized_count, nsamples=nsamples))

if __name__ == "__main__":
    print("i'm deprecated")
    exit(1)
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=CustomFormatter)
    parser.add_argument("-d", "--database", default="waldb6",
                        choices=["waldb", "waldb2", "waldb4", "waldb1", "waldb6"],
                        help="the database to initialize samples in")
    parser.add_argument("--level", default="DEBUG", action=DereferenceKeyAction,
                        choices=LOGGING_LEVELS, help="the logging level to use")
    args = parser.parse_args()
    initialize_samples(args.database, args.level)
