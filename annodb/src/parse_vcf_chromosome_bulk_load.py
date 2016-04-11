#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Take one chromosome from an annotated VCF, create output tables to subsequently
load into the database, and bulk load the novel variants

Brett Copeland <bc2675@cumc.columbia.edu>
4/2016
"""

import argparse
import sys
import MySQLdb
import zlib
from struct import unpack
import match_indels_chromosome
from ConfigParser import SafeConfigParser
from utilities import *
from dragen_globals import *
from db_statements import *

# the output format for each novel variant entry
NOVEL_VARIANT_OUTPUT_FORMAT = (
    "{" + "}\t{".join(
        ["variant_id", "POS", "REF", "ALT", "rs_number", "transcript_stable_id",
         "effect", "HGVS_c", "HGVS_p", "impact", "polyphen_humdiv",
         "polyphen_humvar", "gene", "indel"]) + "}")

def format_NULL_value(value):
    """convert the specified value to \N for NULL where appropriate for
    outputting and subsequent loading
    """
    return value if value else "\\N"

def output_novel_variant_entry(
    novel_fh, variant_id, POS, REF, ALT, rs_number, indel, impact,
    transcript_stable_id="", effect=None, HGVS_c=None, HGVS_p=None,
    polyphen_humdiv=None, polyphen_humvar=None, gene=None):
    """output a specific novel variant entry to novel_fh
    """
    novel_fh.write(NOVEL_VARIANT_OUTPUT_FORMAT.format(
        variant_id=variant_id, POS=POS, REF=REF, ALT=ALT,
        rs_number=format_NULL_value(rs_number),
        transcript_stable_id=transcript_stable_id,
        effect=format_NULL_value(effect), HGVS_c=format_NULL_value(HGVS_c),
        HGVS_p=format_NULL_value(HGVS_p), impact=impact,
        polyphen_humdiv=format_NULL_value(polyphen_humdiv),
        polyphen_humvar=format_NULL_value(polyphen_humvar),
        gene=format_NULL_value(gene), indel=indel) + "\n")

def output_novel_variant(
    novel_fh, variant_id, CHROM, POS, REF, ALT, original_ALT, rs_number, ANNs):
    """output all entries for the novel variant to novel_fh and increment
    variant_id
    """
    rs_number = "" if rs_number == "." else rs_number
    indel = 1 if len(REF) > 1 or len(ALT) > 1 else 0
    anns = []
    for ann in ANNs.split(","):
        alt_allele, values = ann.split("|", 1)
        if alt_allele == original_ALT:
            anns.append(values.split("|"))
    if anns:
        # primary key is POS, variant_id, transcript_stable_id, effect
        # POS, variant_id is invariant so keep track of transcript_id, effect in
        # order to avoid duplicates, i.e. integrity errors
        annotations = {}
        for (effects, impact, gene, gene_id, feature_type, feature_id,
             transcript_biotype, rank_total, HGVS_c, HGVS_p, cDNA_position,
             CDS_position, protein_position, distance, errors) in anns:
            for effect in effects.split("&"):
                annotations_key = (feature_id, effect)
                if annotations_key in annotations:
                    if "N" in REF:
                        # ignore the fact that SnpEff annotates multiple HGVS_c
                        # when there's an N in the reference allele; just take
                        # the first one
                        continue
                    elif effect == "splice_region_variant":
                        # splice_region_variant may be duplicated by SnpEff with
                        # multiple impacts, e.g. HIGH (when joined with other
                        # effects) and LOW - more accurate to take the LOW
                        # impact, so update the annotations (which are sorted by
                        # impact)
                        annotations[annotations_key] = {
                            "impact":impact, "transcript_stable_id":feature_id,
                            "effect":effect, "HGVS_c":HGVS_c, "HGVS_p":HGVS_p,
                            "gene":gene}
                    else:
                        raise ValueError(
                            "error: duplicated ({transcript_stable_id}, {effect}"
                            ") for {CHROM}-{POS}-{REF}-{ALT}".format(
                                transcript_stable_id=feature_id, effect=effect,
                                CHROM=CHROM, POS=POS, REF=REF, ALT=ALT))
                else:
                    annotations[annotations_key] = {
                        "impact":impact, "transcript_stable_id":feature_id,
                        "effect":effect, "HGVS_c":HGVS_c, "HGVS_p":HGVS_p,
                        "gene":gene}
        for annotation_values in annotations.itervalues():
            output_novel_variant_entry(
                novel_fh, variant_id, POS, REF, ALT, rs_number, indel,
                **annotation_values)
    else:
        raise ValueError(
            "error: {CHROM}-{POS}-{REF}-{ALT} has no SnpEff annotation(s)".
            format(CHROM=CHROM, POS=POS, REF=REF, ALT=ALT))

def get_variant_id(novel_fh, cur, CHROM, POS, REF, ALT, rs_number, ANNs):
    """return the variant_id of the given variant and output it to novel_fh
    if it's novel
    """
    # can safely overwrite REF, but need original ALT in order to match up with
    # SnpEff annotations
    REF, alt, offset = simplify_REF_ALT_alleles(REF, ALT)
    POS += offset
    block_id = POS / BLOCK_SIZE
    cur.execute(VARIANT_EXISTS_QUERY.format(
        CHROM=CHROM, POS=POS, REF=REF, ALT=alt))
    row = cur.fetchone()
    if row:
        return row[0], block_id
    else:
        if len(REF) > 1 or len(ALT) > 1:
            # perform indel matching
            matched_indel_id, matched_block_id = match_indels_chromosome(
                cur, CHROM, POS, REF, ALT)
            if matched_indel_id is not None:
                return matched_indels_id, matched_block_id
        global novel_variant_id
        variant_id = novel_variant_id
        output_novel_variant(
            novel_fh, variant_id, CHROM, POS, REF, alt,
            ALT, rs_number, ANNs)
        novel_variant_id += 1
        return variant_id, block_id

def main(input_vcf_fh, output_directory, sample_id, CHROM, log_fh):
    """parse the chromosome VCF, create tables for loading the variants and
    calls, and load the variants
    """
    log_fh.write(log_output(
        "parsing sample_id {sample_id}'s chromosome {CHROM} VCF".format(
            sample_id=sample_id, CHROM=CHROM) + "\n"))
    config_parser = SafeConfigParser()
    config_parser.read(CNF)
    db = MySQLdb.connect(
        read_default_file=config_parser.get("annodb", "DB_CONFIG_FILE"),
        read_default_group="clientdragen01")
    cur = db.cursor()
    # keep track of what the auto-increment value should be for variant_id
    global novel_variant_id
    cur.execute(GET_MAX_VARIANT_ID.format(CHROM=CHROM))
    novel_variant_id = cur.fetchone()[0]
    novel_variants_file = ("{output_directory}/novel_variant_chr{CHROM}.txt".
                           format(output_directory=output_directory, CHROM=CHROM))
    try:
        with open(novel_variants_file, "w") as novel_fh, \
                open("{output_directory}/called_variant_chr{CHROM}.txt".format(
                    output_directory=output_directory, CHROM=CHROM), "w") as calls_fh:
            for line in input_vcf_fh:
                fields = get_vcf_fields_dict(line)
                if fields["CHROM"] != CHROM:
                    raise ValueError(
                        "error: encountered chromosome {chromosome} when "
                        "{CHROM} was expected".format(
                            chromosome=fields["CHROM"], CHROM=CHROM))
                INFO_dict = create_INFO_dict(fields["INFO"])
                INFO_dict["PASS"] = 1 if fields["FILTER"] == "PASS" else 0
                ALT_alleles = fields["ALT"].split(",")
                nalleles = len(ALT_alleles)
                variant_ids = []
                for ALT_allele in ALT_alleles:
                    variant_ids.append(
                        get_variant_id(
                            novel_fh, cur, CHROM, int(fields["POS"]), fields["REF"],
                            ALT_allele, fields["rs_number"], INFO_dict["ANN"]))
                for variant_stat in (
                    "FS", "MQ", "QD", "ReadPosRankSum", "MQRankSum"):
                    if variant_stat not in INFO_dict:
                        # NULL value for loading
                        INFO_dict[variant_stat] = "\\N"
                call_stats = dict(zip(fields["FORMAT"].split(":"),
                                      fields["call"].split(":")))
                call = {"sample_id":sample_id, "GQ":call_stats["GQ"],
                        "QUAL":fields["QUAL"], "DP_pileup":0}
                if nalleles == 1:
                    call["variant_id"] = variant_ids[0][0]
                    call["block_id"] = variant_ids[0][1]
                    GTs = call_stats["GT"].split("/")
                    if len(GTs) == 2:
                        for x, GT in enumerate(GTs):
                            if GT in VALID_GTS:
                                GTs[x] = int(GT)
                            else:
                                raise ValueError(
                                    "error: invalid GT {GT} @ {CHROM}-{POS}-{REF}"
                                    "-{ALT}".format(GT=call_stats["GT"], **fields))
                    else:
                        raise ValueError(
                            "error: invalid GT {GT} @ {CHROM}-{POS}-{REF}"
                            "-{ALT}".format(GT=calls_stats["GT"], **fields))
                    call["GT"] = sum(GTs)
                    AD_REF, AD_ALT = call_stats["AD"].split(",")
                    call["AD_REF"] = AD_REF
                    call["AD_ALT"] = AD_ALT
                    calls_fh.write(VARIANT_CALL_FORMAT.format(
                        **merge_dicts(call, INFO_dict)) + "\n")
                else:
                    if call_stats["GT"] == "1/2":
                        call["GT"] = 1
                    else:
                        raise ValueError(
                            "error: invalid GT {GT} @ {CHROM}-{POS}-{REF}"
                            "-{ALT}".format(GT=call_stats["GT"], **fields))
                    ADs = call_stats["AD"].split(",")
                    call["AD_REF"] = ADs[0]
                    for x, (variant_id, block_id) in enumerate(variant_ids):
                        calls_fh.write(VARIANT_CALL_FORMAT.format(
                            **merge_dicts(
                                call, {"AD_ALT":ADs[1 + x], "variant_id":variant_id,
                                       "block_id":block_id}, INFO_dict)) + "\n")
        log_fh.write(log_output(
            "successfully finished parsing sample_id {sample_id}'s chromosome "
            "{CHROM} VCF".format(
                sample_id=sample_id, CHROM=CHROM) + "\n"))
        cur.execute(LOAD_TABLE.format(
            table_file=novel_variants_file,
            table_name="variant_chr{CHROM}".format(CHROM=CHROM)))
        db.commit()
        log_fh.write(log_output(
            "successfully loaded sample_id {sample_id}'s chromosome {CHROM} "
            "novel variants".format(sample_id=sample_id, CHROM=CHROM) + "\n"))
    finally:
        if db.open:
            db.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=CustomFormatter, description=__doc__)
    parser.add_argument("INPUT_VCF", type=fh, help="the chromosome VCF to parse")
    parser.add_argument("OUTPUT_DIRECTORY", type=make_directory_if_not_exists,
                        help="the directory to output to")
    parser.add_argument("SAMPLE_ID", type=natural_number,
                        help="the sample_id of the sample whose VCF this is")
    parser.add_argument("CHROM", choices=[str(chrom) for chrom in range(1, 23)]
                        + ["X", "Y", "MT"],
                        help="the chromosome that is being processed")
    parser.add_argument("--log", type=argparse.FileType("a"),
                        default=sys.stderr, help="log file to write errors to")
    args = parser.parse_args()
    main(args.INPUT_VCF, args.OUTPUT_DIRECTORY,
         args.SAMPLE_ID, args.CHROM, args.log)
