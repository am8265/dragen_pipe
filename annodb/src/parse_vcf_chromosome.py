#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Take one chromosome from an annotated VCF, insert novel variants,
and create output tables to subsequently load into the database

Brett Copeland <bc2675@cumc.columbia.edu>
3/2016
"""

import argparse
import sys
import MySQLdb
from time import sleep
from ConfigParser import SafeConfigParser
from utilities import *

# the chromosomes to include
CNF = "/nfs/goldstein/software/dragen/dragen.cnf"
VARIANT_EXISTS_QUERY = """
SELECT variant_id
FROM variant_chr{CHROM}
WHERE POS = {POS} AND REF = "{REF}" AND ALT = "{ALT}"
LIMIT 1
"""
VARIANT_INSERT = """
INSERT INTO variant_chr{CHROM} ({variant_id_placeholder}POS, REF, ALT,
    rs_number, transcript_stable_id, effect, HGVS_c, HGVS_p, impact,
    polyphen_humdiv, polyphen_humvar, gene)
VALUE ({variant_id}{POS}, "{REF}", "{ALT}", {rs_number},
    {transcript_stable_id}, {effect}, {HGVS_c}, {HGVS_p}, {impact},
    {polyphen_humdiv}, {polyphen_humvar}, {gene})
"""
VCF_COLUMNS = ["CHROM", "POS", "rs_number", "REF", "ALT", "QUAL",
           "FILTER", "INFO", "FORMAT", "call"]
BLOCK_SIZE = 10000
VARIANT_CALL_FORMAT = (
    "{sample_id}\t{variant_id}\t{block_id}\t{GT}\t{DP_pileup}\t{DP}\t"
    "{AD_REF}\t{AD_ALT}\t{GQ}\t{FS}\t{MQ}\t{QD}\t{QUAL}\t{ReadPosRankSum}\t"
    "{MQRankSum}\t{PASS}")
VALID_GTS = set(["0", "1"])

def get_vcf_fields_dict(line):
    """return a dict of the values of each field in the line
    """
    return dict(zip(VCF_COLUMNS, line.rstrip().split("\t")))

def escaped_value(value):
    """return the properly escaped value for a MySQL insert
    """
    if value:
        return '"{value}"'.format(value=value)
    else:
        return "NULL"

def insert_variant_entry(
    cur, CHROM, POS, REF, ALT, rs_number, variant_id=None,
    transcript_stable_id="", effect="", HGVS_c="", HGVS_p="", impact="",
    polyphen_humdiv="", polyphen_humvar="", gene=""):
    """perform the actual insert and escape any fields appropriately
    """
    insert = (VARIANT_INSERT.format(
        CHROM=CHROM, variant_id_placeholder="variant_id, " if variant_id
        else "", variant_id="{variant_id}, ".format(variant_id=variant_id) if
        variant_id else "", POS=POS, REF=REF, ALT=ALT,
        rs_number=escaped_value(rs_number),
        transcript_stable_id=escaped_value(transcript_stable_id),
        effect=escaped_value(effect),
        HGVS_c=escaped_value(HGVS_c),
        HGVS_p=escaped_value(HGVS_p),
        impact=escaped_value(impact),
        polyphen_humdiv=escaped_value(polyphen_humdiv),
        polyphen_humvar=escaped_value(polyphen_humvar),
        gene=escaped_value(gene)))
    print(insert)
    cur.execute(insert)

def insert_novel_variant(
    db, cur, CHROM, POS, REF, ALT, original_ALT, rs_number, ANNs):
    """insert one or more entries into the variant table for the novel variant
    and return the variant_id
    """
    rs_number = "" if rs_number == "." else rs_number
    anns = []
    for ann in ANNs.split(","):
        alt_allele, values = ann.split("|", 1)
        if alt_allele == original_ALT:
            anns.append(values.split("|"))
    # need to ensure the effects are sorted so that if two jobs try to insert
    # the same variant at the same time they'll process the same effect first
    # and cause a MySQLdb.IntegrityError
    anns = sorted(anns)
    try:
        if anns:
            (effect, impact, gene, gene_id, feature_type, feature_id,
             transcript_biotype, rank_total, HGVS_c, HGVS_p, cDNA_position,
             CDS_position, protein_position, distance, errors) = anns[0]
            effs = effect.split("&")
            insert_variant_entry(
                cur, CHROM, POS, REF, ALT, rs_number,
                transcript_stable_id=feature_id, effect=effs[0], HGVS_c=HGVS_c,
                HGVS_p=HGVS_p, impact=impact, gene=gene)
            cur.execute("SELECT LAST_INSERT_ID()")
            variant_id = cur.fetchone()[0]
            for eff in effs[1:]:
                insert_variant_entry(
                    cur, CHROM, POS, REF, ALT, rs_number, variant_id = variant_id,
                    transcript_stable_id=feature_id, effect=effs[0], HGVS_c=HGVS_c,
                    HGVS_p=HGVS_p, impact=impact, gene=gene)
        else:
            insert_variant_entry(
                cur, CHROM, POS, REF, ALT, rs_number)
            cur.execute("SELECT LAST_INSERT_ID()")
            variant_id = cur.fetchone()[0]
    except MySQLdb.IntegrityError:
        # another sample is currently inserting this variant, so abort trying to
        # insert the variant, wait 10 seconds, and try to get the new variant_id
        sleep(10)
        cur.execute(VARIANT_EXISTS_QUERY.format(
            CHROM=CHROM, POS=POS, REF=REF, ALT=ALT))
        row = cur.fetchone()
        if row:
            # variant is no longer novel, so no need to do the insert
            return row[0]
        else:
            # DB may be unresponsive or there's some other error trying to get
            # this variant_id that should be available
            raise
    # insert the rest of the entries for the variant
    for (effect, impact, gene, gene_id, feature_type, feature_id,
         transcript_biotype, rank_total, HGVS_c, HGVS_p, cDNA_position,
         CDS_position, protein_position, distance, errors) in anns[1:]:
        try:
            # SnpEff concatenates multiple effects with an ampersand; we split
            # these for efficiency querying specific effects
            # alternative is to use InnoDB table with FULLTEXT index
            for eff in effect.split("&"):
                insert_variant_entry(
                    cur, CHROM, POS, REF, ALT, rs_number, variant_id=variant_id,
                    transcript_stable_id=feature_id, effect=eff, HGVS_c=HGVS_c,
                    HGVS_p=HGVS_p, impact=impact, gene=gene)
        except:
            db.rollback()
            raise
    # commit immediately once all annotations are successfully loaded for a
    # variant
    db.commit()
    return variant_id

def get_variant_id(db, cur, CHROM, POS, REF, ALT, rs_number, ANNs):
    """return the variant_id of the given variant and insert if it's novel
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
        return (insert_novel_variant(db, cur, CHROM, POS, REF, alt, ALT,
                                     rs_number, ANNs), block_id)

def main(input_vcf_fh, output_directory, sample_id, CHROM, log_fh):
    """parse the chromosome VCF, insert novel variants to the database,
    and create tables for later loading
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
    try:
        with open("{output_directory}/called_variant_chr{CHROM}.txt".format(
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
                        get_variant_id(db, cur, CHROM, int(fields["POS"]),
                                       fields["REF"], ALT_allele,
                                       fields["rs_number"], INFO_dict["ANN"]))
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
                            "-{ALT}".format(GT=call_stats["GT"], **fields))
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
