#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Truncate all tables for re-importing.
"""

from waldb_globals import get_cfg, get_connection

cfg = get_cfg()
db = get_connection("waldb")
try:
    cur = db.cursor()
    cur.execute("TRUNCATE TABLE matched_indels")
    for CHROM in [
        chromosome[0].upper() for chromosome in cfg.items("chromosomes")]:
        for table_prefix in ("DP_bins_chr", "GQ_bins_chr", "called_variant_chr",
                             "custom_transcript_ids_chr", "variant_chr",
                             "indel_chr"):
            cur.execute("TRUNCATE TABLE {table_prefix}{CHROM}".format(
                table_prefix=table_prefix, CHROM=CHROM))
        cur.execute("ALTER TABLE variant_chr{CHROM} AUTO_INCREMENT=1".format(
            CHROM=CHROM))
    db.commit()
finally:
    if db.open:
        db.close()
