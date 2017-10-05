#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Confirm all knownvar variants are loaded into the database
"""

import argparse
import logging
import match_indels
from collections import Counter
from waldb_globals import *

def main(knownvar_vcf_fh, database_group):
    level = logging.DEBUG
    logger = logging.getLogger(__name__)
    logger.setLevel(level)
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(level)
    formatter = logging.Formatter(cfg.get("logging", "format"))
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    parse_vcf_logger = logging.getLogger("parse_vcf")
    parse_vcf_logger.setLevel(logging.INFO)
    parse_vcf_handler = logging.StreamHandler(sys.stderr)
    parse_vcf_handler.setLevel(logging.INFO)
    parse_vcf_handler.setFormatter(formatter)
    parse_vcf_logger.addHandler(handler)
    db = get_connection(database_group)
    missing_variants = Counter()
    try:
        cur = db.cursor()
        for line in knownvar_vcf_fh:
            if line.startswith("#CHROM"):
                break
        for line in knownvar_vcf_fh:
            fields = line.strip().split("\t")
            chrom, pos = fields[:2]
            ref, alt = fields[3:5]
            ref_sequence = match_indels.get_allele_in_reference_genome(
                chrom, int(pos), ref)
            if ref_sequence != ref:
                logger.debug("Reference mismatch for {chrom}-{pos}-{ref}-{alt}, "
                             "skipping".format(chrom=chrom, pos=pos, ref=ref, alt=alt))
                continue
            cur.execute("SELECT 1 FROM variant_chr{chrom} WHERE POS = {pos} "
                        "AND REF = '{ref}' AND ALT = '{alt}'".format(
                            chrom=chrom, pos=pos, ref=ref, alt=alt))
            row = cur.fetchone()
            if not row:
                cur.execute("SELECT variant_id FROM matched_indels WHERE CHROM = '{chrom}'"
                            " AND POS = {pos} AND REF = '{ref}' AND ALT = '{alt}'".
                            format(chrom=chrom, pos=pos, ref=ref, alt=alt))
                row = cur.fetchone()
                if row:
                    variant_id = row[0]
                    cur.execute("SELECT POS, REF, ALT from variant_chr{chrom} "
                                "WHERE variant_id = {variant_id}".format(
                                    chrom=chrom, variant_id=variant_id))
                    POS, REF, ALT = cur.fetchone()
                    logger.debug("Found matched indel: {chrom}-{pos}-{ref}-{alt} "
                                 "with {chrom}-{POS}-{REF}-{ALT}".format(
                                     chrom=chrom, pos=pos, ref=ref, alt=alt,
                                     POS=POS, REF=REF, ALT=ALT))
                else:
                    logger.warning("Failed to find variant: {chrom}-{pos}-{ref}-{alt}".
                                   format(chrom=chrom, pos=pos, ref=ref, alt=alt))
                    missing_variants[chrom] += 1
    except Exception, m:
        import ipdb
        ipdb.set_trace()
        ipdb.pm()
    finally:
        for chrom in CHROMs:
            logger.info("{} variants missing from chromosome {}".format(
                missing_variants[chrom], chrom))
        if db.open:
            db.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=CustomFormatter, description=__doc__)
    parser.add_argument("KNOWNVAR_VCF", type=argparse.FileType("r"),
                        help="the knownvar VCF to check")
    parser.add_argument("--database", default="waldb6",
                        help="the database to connect to")
    args = parser.parse_args()
    main(args.KNOWNVAR_VCF, args.database)
