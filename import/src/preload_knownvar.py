#!/nfs/goldstein/software/python2.7.7/bin/python

"""
Annotate all Knownvar indels with dbSNP + SnpEff, then perform indel matching +
loading unique indels in order to ensure that all future indels will be in
concordance with the public clinical indel data sets
"""

import argparse
import subprocess
import waldb_globals
import MySQLdb
import logging
import os
import sys
import shlex
import parse_vcf
from collections import defaultdict, OrderedDict

cfg = waldb_globals.get_cfg()

JAVA = "/nfs/goldstein/software/jdk1.8.0_05/bin/java"
CLINEFF = "/nfs/seqscratch09/clineff/clinEff/ClinEff.jar"
DBSNP = "/nfs/seqscratch11/rp2801/ethnicity_1000g/dbsnp_fixed.vcf"
ANNOTATE_VCF = (
    "{java} -Xmx6G -jar {clineff} GRCh37.87 -db {dbsnp} {vcf}")

def main(vcf_fn, level=logging.DEBUG, load_tables=True):
    logger = logging.getLogger(__name__)
    logger.setLevel(level)
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(level)
    formatter = logging.Formatter(cfg.get("logging", "format"))
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    parse_vcf_logger = logging.getLogger("parse_vcf")
    parse_vcf_logger.setLevel(level)
    parse_vcf_handler = logging.StreamHandler(sys.stderr)
    parse_vcf_handler.setLevel(level)
    parse_vcf_handler.setFormatter(formatter)
    parse_vcf_logger.addHandler(handler)
    match_indels_logger = logging.getLogger("match_indels")
    match_indels_logger.setLevel(level)
    match_indels_handler = logging.StreamHandler(sys.stderr)
    match_indels_handler.setLevel(level)
    match_indels_handler.setFormatter(formatter)
    match_indels_logger.addHandler(handler)
    
    output_base = os.path.splitext(vcf_fn)[0] + ".final"
    final_vcf = output_base + ".vcf"
    novel_variants = output_base + ".novel_variants.{CHROM}.txt"
    novel_indels = output_base + ".novel_indels.{CHROM}.txt"
    novel_transcripts = output_base + ".novel_transcripts.{CHROM}.txt"
    matched_indels = output_base + ".matched_indels.{CHROM}.txt"

    clineff_cmd = ANNOTATE_VCF.format(
        java=JAVA, clineff=CLINEFF, dbsnp=DBSNP, vcf=vcf_fn)
    if os.path.exists(final_vcf):
        logger.info("SKIPPING:\n" + clineff_cmd)
    else:
        logger.info(clineff_cmd)
        with open(final_vcf, "w") as out_fh:
            p = subprocess.Popen(shlex.split(clineff_cmd), stdout=out_fh)
            p.communicate()

    db = waldb_globals.get_connection("waldb4")
    cur = db.cursor()
    (effect_rankings, high_impact_effect_ids, moderate_impact_effect_ids,
     low_impact_effect_ids, modifier_impact_effect_ids) = (
         parse_vcf.get_effect_rankings(cur))
    current_chrom = -1
    chrom_files = OrderedDict()

    with open(final_vcf) as vcf_fh:
        for line in vcf_fh:
            if line.startswith("#CHROM"):
                break
        for line in vcf_fh:
            fields = waldb_globals.VCF_fields_dict(line.strip().split("\t"))
            if fields["CHROM"] != current_chrom:
                if current_chrom != -1:
                    novel_variants_fh.close()
                    novel_indels_fh.close()
                    novel_transcripts_fh.close()
                    matched_indels_fh.close()
                novel_variants_fn = novel_variants.format(**fields)
                novel_variants_fh = open(novel_variants_fn, "w")
                novel_indels_fn = novel_indels.format(**fields)
                novel_indels_fh = open(novel_indels_fn, "w")
                novel_transcripts_fn = novel_transcripts.format(**fields)
                novel_transcripts_fh = open(novel_transcripts_fn, "w")
                matched_indels_fn = matched_indels.format(**fields)
                matched_indels_fh = open(matched_indels_fn, "w")
                custom_transcript_ids, novel_transcripts_id = (
                    parse_vcf.get_custom_transcript_ids(cur, fields["CHROM"]))
                current_chrom = fields["CHROM"]
                logger.debug("starting chromosome {}".format(current_chrom))
                chrom_files[current_chrom] = (
                    novel_variants_fn, novel_indels_fn, novel_transcripts_fn,
                    matched_indels_fn)
                novel_variant_id = parse_vcf.get_max_variant_id(
                    cur, current_chrom)

            POS = int(fields["POS"])
            INFO = waldb_globals.create_INFO_dict(fields["INFO"])
            (variant_id, highest_impact, block_id, novel_variant_id,
             novel_transcripts_id) = parse_vcf.get_variant_id(
                 novel_variants_fh, novel_indels_fh, novel_transcripts_fh,
                 matched_indels_fh, cur, 0, current_chrom, POS, fields["REF"],
                 fields["ALT"], fields["rs_number"], INFO["ANN"],
                 novel_variant_id, novel_transcripts_id, effect_rankings,
                 high_impact_effect_ids, moderate_impact_effect_ids,
                 low_impact_effect_ids, modifier_impact_effect_ids, None, None,
                 False, custom_transcript_ids, parse_vcf_logger)
    for fh in (novel_variants_fh, novel_indels_fh, novel_transcripts_fh,
               matched_indels_fh):
        fh.close()

    if load_tables:
        for chrom in chrom_files:
            (novel_variants_fn, novel_indels_fn, novel_transcripts_fn,
             matched_indels_fn) = chrom_files[chrom]
            for table, fn in (("variant_chr{}", novel_variants_fn),
                              ("indel_chr{}", novel_indels_fn),
                              ("custom_transcript_ids_chr{}",
                               novel_transcripts_fn),
                              ("matched_indels", matched_indels_fn)):
                i = "LOAD DATA INFILE '{fn}' INTO TABLE {table}".format(
                    fn=fn, table=table.format(chrom))
                logger.info(i)
                cur.execute(i)
    db.commit()
    db.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=waldb_globals.CustomFormatter)
    parser.add_argument("VCF", type=waldb_globals.file_exists,
                        help="the indel VCF to process")
    parser.add_argument("--level", default="DEBUG",
                        action=waldb_globals.DereferenceKeyAction,
                        choices=waldb_globals.LOGGING_LEVELS,
                        help="the logging level to use")
    parser.add_argument("--dry_run", action="store_true",
                        help="don't actually load the tables generated")
    args = parser.parse_args()
    main(args.VCF, args.level, not args.dry_run)
