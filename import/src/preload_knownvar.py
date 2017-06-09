#!/nfs/goldstein/software/python2.7.7/bin/python

"""
Annotate all Knownvar indels with dbSNP + SnpEff, then perform indel matching +
loading unique indels in order to ensure that all future indels will be in
concordance with the public clinical indel data sets
"""

import argparse
import subprocess
import MySQLdb
import logging
import os
import sys
import shlex
import parse_vcf
from collections import defaultdict, OrderedDict
from waldb_globals import *

cfg = get_cfg()

JAVA = "/nfs/goldstein/software/jdk1.8.0_05/bin/java"
CLINEFF = "/nfs/seqscratch09/clineff/clinEff/ClinEff.jar"
DBSNP = "/nfs/seqscratch11/rp2801/ethnicity_1000g/dbsnp_fixed.vcf"
ANNOTATE_VCF = (
    "{java} -Xmx6G -jar {clineff} GRCh37.87 -db {dbsnp} {vcf}")

def directory(arg):
    d = os.path.realpath(arg)
    if not os.path.isdir(d):
        os.makedirs(d)
    return d

def main(vcf_fn, output_directory, level=logging.DEBUG, load_tables=True):
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
    match_indels_logger = logging.getLogger("match_indels")
    match_indels_logger.setLevel(level)
    match_indels_handler = logging.StreamHandler(sys.stderr)
    match_indels_handler.setLevel(level)
    match_indels_handler.setFormatter(formatter)
    match_indels_logger.addHandler(handler)
    
    output_base = os.path.join(
        output_directory, os.path.splitext(os.path.basename(vcf_fn))[0] + ".final")
    final_vcf = output_base + ".vcf"
    novel_variants = output_base + ".novel_variants.{CHROM}.txt"
    novel_indels = output_base + ".novel_indels.{CHROM}.txt"
    novel_transcripts = output_base + ".novel_transcripts.{CHROM}.txt"
    matched_indels = output_base + ".matched_indels.{CHROM}.txt"
    updated_indels = output_base + ".indels.updated.{CHROM}.txt"

    clineff_cmd = ANNOTATE_VCF.format(
        java=JAVA, clineff=CLINEFF, dbsnp=DBSNP, vcf=vcf_fn)
    if os.path.exists(final_vcf):
        logger.info("SKIPPING:\n" + clineff_cmd)
    else:
        logger.info(clineff_cmd)
        with open(final_vcf, "w") as out_fh:
            p = subprocess.Popen(shlex.split(clineff_cmd), stdout=out_fh)
            p.communicate()

    db = get_connection("waldb4")
    cur = db.cursor()
    (effect_rankings, high_impact_effect_ids, moderate_impact_effect_ids,
     low_impact_effect_ids, modifier_impact_effect_ids) = (
         parse_vcf.get_effect_rankings(cur))
    current_chrom = -1
    chrom_files = OrderedDict()
    polyphen_matrixes_by_stable_id = {"humvar":{}, "humdiv":{}}
    polyphen_stable_ids_to_ignore = {"humvar":set(), "humdiv":set()}
    # delete previous matched_indel records that may have existed if we are
    # updating indel annotations
    matched_indels_to_delete = defaultdict(set)

    with open(final_vcf) as vcf_fh:
        for line in vcf_fh:
            if line.startswith("#CHROM"):
                break
        for line in vcf_fh:
            fields = VCF_fields_dict(line.strip().split("\t"))
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
                updated_indels_fn = updated_indels.format(**fields)
                updated_indels_fh = open(updated_indels_fn, "w")
                custom_transcript_ids, novel_transcripts_id = (
                    parse_vcf.get_custom_transcript_ids(cur, fields["CHROM"]))
                current_chrom = fields["CHROM"]
                logger.debug("starting chromosome {}".format(current_chrom))
                chrom_files[current_chrom] = (
                    novel_variants_fn, novel_indels_fn, novel_transcripts_fn,
                    matched_indels_fn, updated_indels_fn)
                novel_variant_id = parse_vcf.get_max_variant_id(
                    cur, current_chrom)
                indels_already_seen = set()

            POS = int(fields["POS"])
            INFO = create_INFO_dict(fields["INFO"])
            # retain the set of variant_ids that match exactly for indels so we
            # know not to update their metadata if there's a redundancy
            # otherwise, if we observe a matched indel, we'll set the variant
            # annotations to match those from knownvar, put the old
            # position/ref/alt in the matched_indels table, and delete the new
            # position/ref/alf from the matched_indels table if present
            (variant_id, highest_impact, block_id, novel_variant_id,
             novel_transcripts_id) = parse_vcf.get_variant_id(
                 novel_variants_fh, novel_indels_fh, novel_transcripts_fh,
                 matched_indels_fh, cur, 0, current_chrom, POS, fields["REF"],
                 fields["ALT"], fields["rs_number"], INFO["ANN"],
                 novel_variant_id, novel_transcripts_id, effect_rankings,
                 high_impact_effect_ids, moderate_impact_effect_ids,
                 low_impact_effect_ids, modifier_impact_effect_ids,
                 polyphen_matrixes_by_stable_id, polyphen_stable_ids_to_ignore,
                 False, custom_transcript_ids, parse_vcf_logger, True,
                 indels_already_seen)
            if len(fields["REF"]) != len(fields["ALT"]):
                if variant_id not in indels_already_seen:
                    cur.execute("SELECT POS, REF, ALT FROM variant_chr{CHROM} "
                                "WHERE variant_id = {variant_id}".format(
                                    CHROM=current_chrom, variant_id=variant_id))
                    row = cur.fetchone()
                    if row and row[0] != POS:
                        # the indel matched an indel already in the database,
                        # but we want to update annotations to this one, and we
                        # will add the current one to the matched_indels table
                        logger.debug("Updating variant_id {variant_id} from "
                                     "{CHROM}-{POS}-{REF}-{ALT} to "
                                     "{CHROM}-{pos}-{ref}-{alt}".format(
                                         variant_id=variant_id,
                                         CHROM=current_chrom, POS=row[0],
                                         REF=row[1], ALT=row[2], pos=POS,
                                         ref=fields["REF"], alt=fields["ALT"]))
                        matched_indels_fh.write(MATCHED_INDEL_OUTPUT_FORMAT.format(
                            CHROM=current_chrom, variant_id=variant_id,
                            POS=row[0], REF=row[1], ALT=row[2], sample_id=0) + "\n")
                        updated_indels_fh.write(NOVEL_INDEL_OUTPUT_FORMAT.format(
                            variant_id=variant_id, POS=POS,
                            REF=fields["REF"], ALT=fields["ALT"],
                            indel_length=len(fields["ALT"]) - len(fields["REF"]))
                            + "\n")
                        # we will delete a record in matched_indels if it exists
                        # for this tuple and delete all records for this
                        # variant_id
                        matched_indels_to_delete[current_chrom].add(
                            (variant_id, POS))
                    indels_already_seen.add(variant_id)
    for fh in (novel_variants_fh, novel_indels_fh, novel_transcripts_fh,
               matched_indels_fh, updated_indels_fh):
        fh.close()

    if load_tables:
        try:
            for chrom in chrom_files:
                (novel_variants_fn, novel_indels_fn, novel_transcripts_fn,
                 matched_indels_fn, updated_indels_fn) = chrom_files[chrom]
                nvariants_to_delete = len(matched_indels_to_delete[chrom])
                logger.info("Deleting {nvariants} indels from variant, indel, and "
                            "matched_indels tables for chromosome {chrom}".format(
                                nvariants=nvariants_to_delete, chrom=chrom))
                for variant_id, POS in matched_indels_to_delete[chrom]:
                    logger.debug("Deleting variant_id {variant_id} on chromosome {chrom}".
                                 format(variant_id=variant_id, chrom=chrom))
                    for table in ("variant_chr", "indel_chr"):
                        cur.execute("DELETE FROM {table}{chrom} WHERE variant_id "
                                    "= {variant_id}".format(
                                        table=table, chrom=chrom, variant_id=variant_id))
                    cur.execute("DELETE FROM matched_indels WHERE CHROM = '{chrom}' "
                                "AND variant_id = {variant_id} AND POS = {POS}".format(
                                    chrom=chrom, variant_id=variant_id, POS=POS))
                for table, fn in (("variant_chr{}", novel_variants_fn),
                                  ("indel_chr{}", novel_indels_fn),
                                  ("custom_transcript_ids_chr{}",
                                   novel_transcripts_fn),
                                  ("matched_indels", matched_indels_fn),
                                   ("indel_chr{}", updated_indels_fn)):
                    i = "LOAD DATA INFILE '{fn}' INTO TABLE {table}".format(
                        fn=fn, table=table.format(chrom))
                    logger.info(i)
                    cur.execute(i)
                logger.debug("Confirming presence of new variants in chromosome "
                             "{chrom}".format(chrom=chrom))
                for variant_id, _ in matched_indels_to_delete[chrom]:
                    for table in ("variant_chr", "indel_chr"):
                        cur.execute("SELECT 1 FROM {table}{chrom} WHERE variant_id "
                                    "= {variant_id}".format(
                                        table=table, chrom=chrom,
                                        variant_id=variant_id))
                        row = cur.fetchone()
                        if not row:
                            # we must have deleted a record from the variant/indel
                            # table and not replaced it; better roll back and debug
                            db.rollback()
                            raise MySQLdb.IntegrityError(
                                "Missing variant_id {variant_id} in {table}{chrom}".
                                format(variant_id=variant_id, table=table, chrom=chrom))
        except Exception, m:
            logger.warning(m)
            db.rollback()
    db.commit()
    db.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=CustomFormatter)
    parser.add_argument("VCF", type=file_exists,
                        help="the indel VCF to process")
    parser.add_argument("-o", "--output_directory", type=directory,
                        default=os.getcwd(), help="the output directory (make "
                        "sure to output to a directory that the DB can load "
                        "from if actually loading data)")
    parser.add_argument("--level", default="DEBUG",
                        action=DereferenceKeyAction,
                        choices=LOGGING_LEVELS,
                        help="the logging level to use")
    parser.add_argument("--dry_run", action="store_true",
                        help="don't actually load the tables generated")
    args = parser.parse_args()
    main(args.VCF, args.output_directory, args.level, not args.dry_run)
