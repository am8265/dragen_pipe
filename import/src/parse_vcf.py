#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Parse a single sample's VCF for a single chromosome and
create output tables to load into the database
"""

import zlib
from struct import unpack
import match_indels_chromosome
import re
from dragen_globals import *
from db_statements import *
from collections import defaultdict
from functools import partial
from time import time
import tabix
import logging

cfg = get_cfg()
HGVS_P_REGEX = re.compile(HGVS_P_PATTERN)

def format_NULL_value(value):
    """convert the specified value to \N for NULL where appropriate for
    outputting and subsequent loading
    """
    return value if value else "\\N"

def call_is_high_quality(QUAL, MQ, FILTER, DP=None):
    """return whether the given call is high quality
    """
    high_quality = (QUAL >= 30 and MQ >= 40
            and FILTER in ("PASS", "LIKELY", "INTERMEDIATE"))
    if high_quality and DP is not None:
        # there is no DP if MNP
        high_quality = high_quality and DP >= 3
    return high_quality

def get_custom_transcript_ids(cur, CHROM):
    """return a dict of all novel transcript ids and the new value
    to begin using
    """
    custom_transcript_ids = {}
    cur.execute(GET_CUSTOM_TRANSCRIPT_IDS.format(CHROM=CHROM))
    for transcript_id, transcript_ids in cur.fetchall():
        custom_transcript_ids[transcript_ids] = transcript_id
    cur.execute(GET_MIN_CUSTOM_TRANSCRIPT_ID.format(CHROM=CHROM))
    min_transcript_id = cur.fetchone()[0]
    if min_transcript_id is None:
        min_transcript_id = 0
    return (custom_transcript_ids, min_transcript_id - 1)

def calculate_polyphen_scores(
    cur, transcript_stable_id, HGVS_p, VariantID,
    polyphen_matrixes_by_stable_id, polyphen_stable_ids_to_ignore, logger):
    """return the PolyPhen scores for the given missense variant
    """
    # cache the transcripts' matrixes for those that have them and cache the
    # transcript IDs for those that don't to avoid having to re-query
    hgvs_p_match = HGVS_P_REGEX.match(HGVS_p)
    scores = {}
    if hgvs_p_match:
        d = hgvs_p_match.groupdict()
        codon_position = int(d["codon_position"])
        amino_acid_changes = d["amino_acid_changes"]
        len_amino_acid_changes = len(amino_acid_changes)
        if len_amino_acid_changes % 3:
            raise ValueError("{VariantID} has an amino acid change length not "
                             "divisible by 3: {HGVS_p}".format(
                                 VariantID=VariantID, HGVS_p=HGVS_p))
        num_amino_acid_changes = len_amino_acid_changes / 3
        if num_amino_acid_changes != 0:
            logger.debug("Encountered multiple impacted amino acids for "
                         "{VariantID}: {HGVS_p}".format(
                             VariantID=VariantID, HGVS_p=HGVS_p))
        max_polyphen_scores = {"humdiv":None, "humvar":None}
        for offset in xrange(num_amino_acid_changes):
            matrix_offset = (3 + 2 * ((codon_position + offset - 1) * 20 +
                                      AMINO_ACIDS[amino_acid_changes[
                                          offset * 3:(offset + 1) * 3]]))
            for polyphen_score in ("humdiv", "humvar"):
                if (transcript_stable_id in
                    polyphen_stable_ids_to_ignore[polyphen_score]):
                    scores["polyphen_" + polyphen_score] = None
                    continue
                if (transcript_stable_id not in
                    polyphen_matrixes_by_stable_id[polyphen_score]):
                    cur.execute(GET_TRANSLATION_MD5_ID.format(
                        stable_id=transcript_stable_id))
                    md5_id_row = cur.fetchone()
                    if md5_id_row:
                        translation_md5_id = md5_id_row[0]
                    else:
                        logger.debug("Could not find {transcript_stable_id}'s "
                                     "MD5 translation ID".format(
                                         transcript_stable_id=transcript_stable_id))
                        polyphen_stable_ids_to_ignore[polyphen_score].add(
                            transcript_stable_id)
                        scores["polyphen_" + polyphen_score] = None
                        continue
                    cur.execute(GET_POLYPHEN_PREDICTION_MATRIX.format(
                        translation_md5_id=translation_md5_id,
                        attrib_id=POLYPHEN_ATTRIB_ID[polyphen_score]))
                    polyphen_matrix_row = cur.fetchone()
                    if polyphen_matrix_row:
                        polyphen_matrixes_by_stable_id[polyphen_score][
                            transcript_stable_id] = (
                                zlib.decompress(polyphen_matrix_row[0],
                                                16 + zlib.MAX_WBITS))
                    else:
                        logger.debug("Could not find {transcript_stable_id}'s "
                                     "PolyPhen-2 {score} matrix".format(
                                         transcript_stable_id=transcript_stable_id,
                                         score=polyphen_score))
                        polyphen_stable_ids_to_ignore[polyphen_score].add(
                            transcript_stable_id)
                        scores["polyphen_" + polyphen_score] = None
                        continue
                    packed_score = polyphen_matrixes_by_stable_id[polyphen_score][
                        transcript_stable_id][matrix_offset:matrix_offset + 2]
                    unpacked_value = unpack("H", packed_score)[0]
                    prediction = int(unpacked_value >> 14)
                    if prediction == 3:
                        # encodes "UNKNOWN" score, so we'll just store NULL
                        continue
                    max_polyphen_scores[polyphen_score] = max(
                        max_polyphen_scores[polyphen_score],
                        unpacked_value & POLYPHEN_PROB_BITMASK)
        # set to the max in the affected interval to accommodate for an MNP
        # potentially altering multiple amino acids
        scores["polyphen_humdiv"] = max_polyphen_scores["humdiv"]
        scores["polyphen_humvar"] = max_polyphen_scores["humvar"]
    else:
        raise ValueError(
            "error: could not parse HGVS_p {HGVS_p} for "
            "{VariantID}".format(
                HGVS_p=HGVS_p, VariantID=VariantID))
    return scores

def get_variant_id(novel_fh, novel_transcripts_fh, matched_indels_fh, cur,
                   CHROM, POS, REF, ALT, rs_number, ANNs, novel_variant_id,
                   novel_transcripts_id, effect_rankings,
                   high_impact_effect_ids, moderate_impact_effect_ids,
                   low_impact_effect_ids, modifier_impact_effect_ids,
                   polyphen_matrixes_by_stable_id,
                   polyphen_stable_ids_to_ignore, high_quality_call,
                   custom_transcript_ids, logger):
    """return the variant_id of the given variant and output it to novel_fh
    if it's novel
    """
    # can safely overwrite REF, but need original ALT in order to match up with
    # SnpEff annotations
    REF, alt, offset = simplify_REF_ALT_alleles(REF, ALT)
    indel_length = len(alt) - len(REF)
    if len(REF) > 255:
        REF = REF[:255]
    if len(alt) > 255:
        alt = alt[:255]
    POS += offset
    block_id = POS / BLOCK_SIZE
    cur.execute(VARIANT_EXISTS_QUERY.format(
        CHROM=CHROM, POS=POS, REF=REF, ALT=alt, indel_length=indel_length))
    rows = cur.fetchall()
    if rows:
        variant_id = rows[0][0]
        effect_ids = [row[1] for row in rows]
        has_high_quality_call =  rows[0][2]
        # treat the variant as novel if it doesn't have a high quality call in
        # the DB and the new call is high quality, so as to update the field
        novel = not has_high_quality_call and high_quality_call
        update_novel_variant_id = False
    else:
        novel = True
        update_novel_variant_id = True
        if indel_length:
            # don't treat as an indel if the length of both is the same, i.e.
            # it's an MNV
            # perform indel matching
            matched_indel_id, matched_block_id = match_indels_chromosome.match_indel(
                cur, CHROM, POS, REF, alt, indel_length)
            if matched_indel_id is not None:
                variant_id = matched_indel_id
                block_id = matched_block_id
                novel = False
                cur.execute(GET_VARIANT_EFFECTS.format(
                    CHROM=CHROM, variant_id=variant_id))
                effect_ids = [row[0] for row in cur.fetchall()]
                # check if this matched indel is already in the matched_indels
                # table
                cur.execute(MATCHED_INDEL_EXISTS.format(
                    CHROM=CHROM, POS=POS, REF=REF, ALT=ALT))
                if not cur.fetchone():
                    matched_indels_fh.write(MATCHED_INDEL_OUTPUT_FORMAT.format(
                        CHROM=CHROM, variant_id=variant_id,
                        POS=POS, REF=REF, ALT=ALT))
    if novel:
        if update_novel_variant_id:
            variant_id = novel_variant_id
            novel_variant_id += 1
            if indel_length:
                # add this indel to the matched indel set so in the corner case
                # of a single sample having > 1 redundant indel, we only create
                # one single novel variant
                match_indels_chromosome.add_new_indel(
                    variant_id, CHROM, POS, REF, ALT, indel_length)
        effect_ids, novel_transcripts_id = output_novel_variant(
            novel_fh, novel_transcripts_fh, cur, variant_id,
            novel_transcripts_id, CHROM, POS,
            REF, alt, indel_length, ALT, rs_number, ANNs, effect_rankings,
            polyphen_matrixes_by_stable_id, polyphen_stable_ids_to_ignore,
            high_quality_call, custom_transcript_ids, logger)
    if any(effect_id in high_impact_effect_ids for effect_id in effect_ids):
        highest_impact = "HIGH"
    elif any(effect_id in moderate_impact_effect_ids
             for effect_id in effect_ids):
        highest_impact = "MODERATE"
    elif any(effect_id in low_impact_effect_ids for effect_id in effect_ids):
        highest_impact = "LOW"
    elif any(effect_id in modifier_impact_effect_ids for effect_id in
             effect_ids):
        highest_impact = "MODIFIER"
    else:
        highest_impact = "\\N"
    return (variant_id, highest_impact, block_id,
            novel_variant_id, novel_transcripts_id)

def output_novel_variant(
    novel_fh, novel_transcripts_fh, cur, variant_id, novel_transcripts_id, CHROM,
    POS, REF, ALT, indel_length, original_ALT, rs_number, ANNs, effect_rankings,
    polyphen_matrixes_by_stable_id, polyphen_stable_ids_to_ignore,
    high_quality_call, custom_transcript_ids, logger,
    impact_ordering=["HIGH", "MODERATE", "LOW", "MODIFIER"]):
    """output all entries for the novel variant to novel_fh and increment
    variant_id
    """
    transcript_ids_dict = {}
    VariantID = "{CHROM}-{POS}-{REF}-{ALT}".format(
        CHROM=CHROM, POS=POS, REF=REF, ALT=ALT)
    rs_number = ("" if rs_number == "." else
                 int(strip_prefix(re.split(";|,", rs_number)[0], "rs")))
    indel = 1 if indel_length else 0
    anns = []
    for ann in ANNs.split(","):
        alt_allele, values = ann.split("|", 1)
        if alt_allele == original_ALT:
            anns.append(values.split("|"))
    effect_ids = []
    if anns:
        # primary key is POS, variant_id, transcript_stable_id, effect
        # POS, variant_id is invariant so keep track of transcript_id, effect in
        # order to avoid duplicates, i.e. integrity errors
        annotations = {}
        # keep track of the effects by transcript in order to annotate the
        # custom annotations missense_variant+splice_region_variant and
        # splice_region_variant+synonymous_variant
        annotations_by_transcript = defaultdict(set)
        for (effects, impact, gene, gene_id, feature_type, feature_id,
             transcript_biotype, rank_total, HGVS_c, HGVS_p, cDNA_position,
             CDS_position, protein_position, distance, errors) in anns:
            # there appears to be a bug such that
            # 3_prime_UTR_truncation+exon_loss appears as
            # 3_prime_UTR_truncation&exon_loss and
            # 5_prime_UTR_truncation+exon_loss_variant appears as
            # 5_prime_UTR_truncation&exon_loss_variant
            if "prime_UTR_truncation" in effects:
                effects = (effects.replace("3_prime_UTR_truncation&exon_loss",
                                           "3_prime_UTR_truncation+exon_loss").
                           replace("5_prime_UTR_truncation&exon_loss_variant",
                                   "5_prime_UTR_truncation+exon_loss_variant"))
            if "non_canonical_start_codon" in effects:
                effects = effects.replace(
                    "initiator_codon_variant&non_canonical_start_codon",
                    "initiator_codon_variant+non_canonical_start_codon")
            for effect in effects.split("&"):
                if effect == "custom":
                    # these correspond to the deprecated INTRON_EXON_BOUNDARY
                    # annotations which SnpEff now natively annotates
                    continue
                if effect not in effect_rankings:
                    raise ValueError(
                        "error: invalid effect {effect} for {VariantID}".format(
                        effect=effect, VariantID=VariantID))
                if impact not in effect_rankings[effect]:
                    # this is likely due to SnpEff concatenating two or more
                    # effects into one annotation - we will try to select the
                    # next least deleterious impact that is valid for the effect
                    impact_idx = impact_ordering.index(impact)
                    corrected_impact_found = False
                    for next_impact in impact_ordering[impact_idx + 1:]:
                        if next_impact in effect_rankings[effect]:
                            impact = next_impact
                            corrected_impact_found = True
                            break
                    if not corrected_impact_found:
                        raise ValueError(
                            "failed to find a valid (impact, effect) to match "
                            "({impact}, {effect})".format(
                                impact=impact, effect=effect))
                # sometimes SnpEff can annotate the same effect in transcripts
                # and treat different case differently, but this will cause
                # integrity errors, so they're converted to upper-case always
                feature_id = feature_id.upper()
                if feature_id == "":
                    transcript_ids_dict[""] = 0
                elif "-" in feature_id or feature_id.startswith("ENSG"):
                    if feature_id not in custom_transcript_ids:
                        # new transcript id
                        custom_transcript_ids[feature_id] = novel_transcripts_id
                        novel_transcripts_fh.write("{id}\t{feature_id}\n".format(
                            id=novel_transcripts_id, feature_id=feature_id))
                        novel_transcripts_id -= 1
                    transcript_ids_dict[feature_id] = (
                        custom_transcript_ids[feature_id])
                else:
                    if feature_id.startswith("ENST"):
                        transcript_ids_dict[feature_id] = int(feature_id[4:])
                    else:
                        raise ValueError(
                            "failure getting a transcript ID for {}".format(
                                feature_id))
                annotations_key = (feature_id, effect)
                effect_id = effect_rankings[effect][impact]
                effect_ids.append(effect_id)
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
                            "transcript_stable_id":transcript_ids_dict[feature_id],
                            "effect_id":effect_id, "HGVS_c":HGVS_c,
                            "HGVS_p":HGVS_p, "gene":gene}
                    elif effect == "exon_intron_boundary":
                        # this can be repeated for unkown reason(s)
                        # will just skip any after the first
                        continue
                else:
                    annotations[annotations_key] = {
                        "transcript_stable_id":transcript_ids_dict[feature_id],
                        "effect_id":effect_id, "HGVS_c":HGVS_c,
                        "HGVS_p":HGVS_p, "gene":gene}
                    annotations_by_transcript[feature_id].add(effect)
                    if effect == "missense_variant":
                        # calculate PolyPhen scores if possible
                        annotations[annotations_key].update(
                            calculate_polyphen_scores(
                                cur, feature_id, HGVS_p, VariantID,
                                polyphen_matrixes_by_stable_id,
                                polyphen_stable_ids_to_ignore, logger))
        for feature_id, effects in annotations_by_transcript.iteritems():
            if "splice_region_variant" in effects:
                if "missense_variant" in effects:
                    # if a single transcript has both these effects, add this
                    # custom effect annotation
                    effect = "missense_variant+splice_region_variant"
                    annotations[(feature_id, effect)] = (
                        annotations[(feature_id, "missense_variant")]).copy()
                    # need to replace the effect_id for this so as to indicate
                    # the proper effect type
                    effect_id = effect_rankings[effect]["MODERATE"]
                    annotations[(feature_id, effect)]["effect_id"] = effect_id
                    effect_ids.append(effect_id)
                if "synonymous_variant" in effects:
                    effect = "splice_region_variant+synonymous_variant"
                    annotations[(feature_id, effect)] = (
                            annotations[(feature_id, "synonymous_variant")]).copy()
                    effect_id = effect_rankings[effect]["LOW"]
                    annotations[(feature_id, effect)]["effect_id"] = effect_id
                    effect_ids.append(effect_id)
        for annotation_values in annotations.itervalues():
            output_novel_variant_entry(
                novel_fh, variant_id, POS, REF, ALT, rs_number, indel,
                indel_length, high_quality_call, logger, **annotation_values)
    else:
        raise ValueError(
            "error: {VariantID} has no SnpEff annotation(s)".
            format(VariantID=VariantID))
    return effect_ids, novel_transcripts_id

def output_novel_variant_entry(
    novel_fh, variant_id, POS, REF, ALT, rs_number, indel, indel_length,
    high_quality_call, logger, transcript_stable_id="", effect_id=None,
    HGVS_c=None, HGVS_p=None, polyphen_humdiv=None,
    polyphen_humvar=None, gene=None):
    """output a specific novel variant entry to novel_fh
    """
    novel_fh.write(NOVEL_VARIANT_OUTPUT_FORMAT.format(
        variant_id=variant_id, POS=POS, REF=REF, ALT=ALT,
        rs_number=format_NULL_value(rs_number),
        transcript_stable_id=transcript_stable_id,
        effect_id=format_NULL_value(effect_id),
        HGVS_c=format_NULL_value(HGVS_c)[:255],
        HGVS_p=format_NULL_value(HGVS_p)[:255],
        polyphen_humdiv=format_NULL_value(polyphen_humdiv),
        polyphen_humvar=format_NULL_value(polyphen_humvar),
        gene=format_NULL_value(gene), indel=indel,
        indel_length=indel_length,
        has_high_quality_call=int(high_quality_call)) + "\n")

def parse_vcf(vcf, CHROM, sample_id, output_base, chromosome_length=None):
    logger = logging.getLogger(__name__)
    logger.info("starting CHROM {}\n".format(CHROM))
    line_count = 0
    start_time = time()
    cfg = get_cfg()
    # store the set of indel_ids for calls, because in rare circumstances there
    # are multiple entries for the same indel
    # (because of repetitive genome/indel matching)
    indel_ids = set()
    db = get_connection("dragen")
    try:
        cur = db.cursor()
        cur.execute(GET_MAX_VARIANT_ID.format(CHROM=CHROM))
        max_variant_id = cur.fetchone()[0]
        novel_variant_id = max_variant_id
        effect_rankings = defaultdict(lambda: defaultdict(int))
        high_impact_effect_ids = set()
        moderate_impact_effect_ids = set()
        low_impact_effect_ids = set()
        modifier_impact_effect_ids = set()
        cur.execute(GET_EFFECT_RANKINGS)
        for effect_id, impact, effect in cur.fetchall():
            effect_rankings[effect][impact] = effect_id
            if impact == "HIGH":
                high_impact_effect_ids.add(effect_id)
            elif impact == "MODERATE":
                moderate_impact_effect_ids.add(effect_id)
            elif impact == "LOW":
                low_impact_effect_ids.add(effect_id)
            elif impact == "MODIFIER":
                modifier_impact_effect_ids.add(effect_id)
        polyphen_matrixes_by_stable_id = {"humvar":{}, "humdiv":{}}
        polyphen_stable_ids_to_ignore = {"humvar":set(), "humdiv":set()}
        novel_variants = output_base + ".novel_variants.txt"
        novel_transcripts = output_base + ".novel_transcripts.txt"
        calls = output_base + ".calls.txt"
        variant_id_vcf = output_base + ".variant_id.vcf"
        matched_indels = output_base + ".matched_indels.txt"
        custom_transcript_ids, novel_transcripts_id = get_custom_transcript_ids(
            cur, CHROM)
        vcf_tabix = tabix.open(vcf)
        with open(novel_variants, "w") as novel_fh, \
                open(novel_transcripts, "w") as novel_transcripts_fh, \
                open(calls, "w") as calls_fh, \
                open(variant_id_vcf, "w") as vcf_out, \
                open(matched_indels, "w") as matched_indels_fh:
            for x, line_fields in enumerate(vcf_tabix.querys(CHROM)):
                fields = VCF_fields_dict(line_fields)
                if fields["CHROM"] != CHROM:
                    raise ValueError(
                        "error: encountered chromosome {chromosome} when "
                        "{CHROM} was expected".format(
                            CHROM=CHROM, chromosome=fields["CHROM"]))
                POS = int(fields["POS"])
                INFO = create_INFO_dict(fields["INFO"])
                if fields["FILTER"] == "PASS":
                    INFO["FILTER"] = "PASS"
                elif fields["FILTER"] == "INDEL_filter":
                    INFO["FILTER"] = "FAIL"
                elif fields["FILTER"] == "VQSRTrancheSNP90.00to99.00":
                    INFO["FILTER"] = "LIKELY"
                elif fields["FILTER"] == "VQSRTrancheSNP99.00to99.90":
                    INFO["FILTER"] = "INTERMEDIATE"
                elif fields["FILTER"] == "VQSRTrancheSNP99.90to100.00":
                    INFO["FILTER"] = "FAIL"
                elif fields["FILTER"] == "VQSRTrancheINDEL90.00to99.00":
                    INFO["FILTER"] = "LIKELY"
                elif fields["FILTER"] == "VQSRTrancheINDEL99.00to99.90":
                    INFO["FILTER"] = "INTERMEDIATE"
                elif fields["FILTER"] == "VQSRTrancheINDEL99.90to100.00":
                    INFO["FILTER"] = "FAIL"
                else:
                    raise ValueError("invalid FILTER {} @ line {}".format(
                        fields["FILTER"], x))
                call_stats = create_call_dict(fields["FORMAT"], fields["call"])
                call = {"sample_id":sample_id, "GQ":call_stats["GQ"],
                        "QUAL":int(round(float(fields["QUAL"])))}
                if "DP" in call_stats:
                    call["DP"] = call_stats["DP"]
                else:
                    # this should always be present except at MNP sites (for
                    # whatever reason)
                    call["DP"] = "\\N"
                ALT_alleles = fields["ALT"].split(",")
                nalleles = len(ALT_alleles)
                if nalleles > 2:
                    raise ValueError("{CHROM}-{POS} has {nalleles} alternate "
                                     "alleles".format(
                                         CHROM=CHROM, POS=POS, nalleles=len(ALT_alleles)))
                # we truncate these to 255 characters, so if they differ after
                # that we'll get different variants that have the same
                # "alternate" - our solution is to just retain the first entry
                # in such a rare instance
                elif nalleles == 2:
                    if ALT_alleles[0][:255] == ALT_alleles[1][:255]:
                        logger.debug(
                            "{CHROM}-{POS} has two alternate alleles with > 255 "
                            "length that are identical to that point".format(
                                CHROM=CHROM, POS=POS))
                        del ALT_alleles[1]
                        nalleles = 1
                        # need to change the GT to 1/1 instead of 1/2
                        call_stats["GT"] = call_stats["GT"].replace("2", "1")
                        # also need to change AD from x,y,z to x,y+z
                        ads = call_stats["AD"].split(",")
                        ads[1] = str(int(ads[1]) + int(ads[2]))
                        del ads[2]
                        call_stats["AD"] = ",".join(ads)
                high_quality_call = call_is_high_quality(
                    int(round(float(fields["QUAL"]))),
                    int(round(float(INFO["MQ"]))) if "MQ" in INFO else 0,
                    INFO["FILTER"],
                    DP=int(call_stats["DP"]) if "DP" in call_stats else None)
                variant_ids = []
                for ALT_allele in ALT_alleles:
                    (variant_id, highest_impact, block_id, novel_variant_id,
                     novel_transcripts_id) = get_variant_id(
                         novel_fh, novel_transcripts_fh, matched_indels_fh,
                         cur, CHROM, POS, fields["REF"], ALT_allele,
                         fields["rs_number"], INFO["ANN"], novel_variant_id,
                         novel_transcripts_id, effect_rankings, high_impact_effect_ids,
                         moderate_impact_effect_ids, low_impact_effect_ids,
                         modifier_impact_effect_ids, polyphen_matrixes_by_stable_id,
                         polyphen_stable_ids_to_ignore, high_quality_call,
                         custom_transcript_ids, logger)
                    indel = len(fields["REF"]) != len(ALT_allele)
                    variant_ids.append(
                        (variant_id, block_id, highest_impact, indel))
                for variant_stat in ("MQ", "QD"):
                    # we store these two values as ints instead of float
                    if variant_stat in INFO:
                        INFO[variant_stat] = int(round(float(INFO[variant_stat])))
                for variant_stat in (
                    "FS", "MQ", "QD", "ReadPosRankSum", "MQRankSum", "VQSLOD"):
                    if variant_stat not in INFO:
                        # NULL value for loading
                        INFO[variant_stat] = "\\N"
                if nalleles == 1:
                    (call["variant_id"], call["block_id"],
                     call["highest_impact"], indel) = variant_ids[0]
                    output_call = True
                    if indel:
                        if call["variant_id"] in indel_ids:
                            # skip outputting this call, because this is a
                            # matched indel which we already output
                            output_call = False
                        else:
                            indel_ids.add(call["variant_id"])
                if nalleles == 1 and output_call:
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
                    if "AD" in call_stats:
                        AD_REF, AD_ALT = call_stats["AD"].split(",")
                        call["AD_REF"] = AD_REF
                        call["AD_ALT"] = AD_ALT
                    else:
                        # these are not present in MNPs
                        call["AD_REF"] = "\\N"
                        call["AD_ALT"] = "\\N"
                    try:
                        gg = VARIANT_CALL_FORMAT.format(
                            **merge_dicts(call, INFO)) + "\n"
                    except KeyError:
                        from pprint import pprint
                        pprint(call)
                        pprint(INFO)
                        raise
                    calls_fh.write(gg)
                elif nalleles != 1:
                    if call_stats["GT"] == "1/2":
                        call["GT"] = 1
                    else:
                        raise ValueError(
                            "error: invalid GT {GT} @ {CHROM}-{POS}-{REF}"
                            "-{ALT}".format(GT=call_stats["GT"], **fields))
                    ADs = call_stats["AD"].split(",")
                    call["AD_REF"] = ADs[0]
                    for x, (variant_id, block_id, highest_impact, indel) in enumerate(variant_ids):
                        if indel:
                            if variant_id in indel_ids:
                                continue
                            else:
                                indel_ids.add(variant_id)
                        calls_fh.write(VARIANT_CALL_FORMAT.format(
                            **merge_dicts(
                                call, {"AD_ALT":ADs[1 + x], "variant_id":variant_id,
                                       "block_id":block_id,
                                       "highest_impact":highest_impact}, INFO)) + "\n")
                # annotate each variant with its variant_id for storage purposes
                line_fields[VCF_COLUMNS_DICT["INFO"]] = (
                    "VariantID=" + ",".join(
                        str(variant_id[0]) for variant_id in variant_ids) +
                    ";" + fields["INFO"])
                vcf_out.write("\t".join(line_fields) + "\n")
        # confirm the highest variant_id is still the same; otherwise something
        # else (e.g. another job added a variant and we can run into collisions)
        # we won't get new results unless we commit/reconnect
        db.commit()
        cur.execute(GET_MAX_VARIANT_ID.format(CHROM=CHROM))
        current_max_variant_id = cur.fetchone()[0]
        if max_variant_id != current_max_variant_id:
            raise ValueError(
                "The max variant_id while starting, {max_variant_id}, does not "
                "match the current max variant_id, {current_max_variant_id}".
                format(max_variant_id=max_variant_id,
                       current_max_variant_id=current_max_variant_id))
    finally:
        logger.info("elapsed time: {} seconds\n".format(
            time() - start_time))
        if db.open:
            db.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=CustomFormatter, description=__doc__)
    parser.add_argument("VCF", type=file_exists,
                        help="the chromosome VCF to parse")
    parser.add_argument("CHROMOSOME",
                        choices=[chrom[0].upper() for chrom in cfg.items("chromosomes")],
                        help="the chromosome that is being processed")
    parser.add_argument("SAMPLE_ID",
                        type=partial(valid_numerical_argument, arg_name="sample_id"),
                        help="the id of the sample")
    parser.add_argument("OUTPUT_BASE", help="the base output file name structure")
    parser.add_argument("-l", "--level", default="ERROR",
                        choices=LOGGING_LEVELS.keys(),
                        help="specify the logging level to use")
    args = parser.parse_args()
    output_base_rp = os.path.realpath(args.OUTPUT_BASE)
    if not os.path.isdir(os.path.dirname(output_base_rp)):
        os.makedirs(os.path.dirname(output_base_rp))
    level = LOGGING_LEVELS[args.level]
    import sys
    logger = logging.getLogger(__name__)
    logger.setLevel(level)
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(level)
    formatter = logging.Formatter(cfg.get("logging", "format"))
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    parse_vcf(args.VCF, args.CHROMOSOME, args.SAMPLE_ID, output_base_rp)
