"""
MySQL statements used in the pipeline
"""

# check novelty of variant/get its id
VARIANT_EXISTS_QUERY = """
SELECT variant_id
FROM variant_chr{CHROM}
WHERE POS = {POS} AND REF = "{REF}" AND ALT = "{ALT}"
LIMIT 1
"""
# insert a novel variant
VARIANT_INSERT = """
INSERT INTO variant_chr{CHROM} ({variant_id_placeholder}POS, REF, ALT,
    rs_number, transcript_stable_id, effect, HGVS_c, HGVS_p, impact,
    polyphen_humdiv, polyphen_humvar, gene, indel)
VALUE ({variant_id}{POS}, "{REF}", "{ALT}", {rs_number},
    {transcript_stable_id}, {effect}, {HGVS_c}, {HGVS_p}, {impact},
    {polyphen_humdiv}, {polyphen_humvar}, {gene}, {indel})
"""
# handle SnpEff "bug" of sometimes outputting multiple impacts for
# splice_region_variants, and take the least impactful (last)
UPDATE_SPLICE_REGION_VARIANT = """
UPDATE variant_chr{CHROM}
SET impact = "{impact}"
WHERE POS = {POS} AND variant_id = {variant_id} AND
    transcript_stable_id = "{transcript_stable_id}" AND effect = "effect"
"""
# get the current maximum variant id in the specified table
GET_MAX_VARIANT_ID = """
SELECT auto_increment
FROM information_schema.tables
WHERE table_name = "variant_chr{CHROM}" AND table_schema = DATABASE()
"""
LOAD_TABLE = """
LOAD DATA INFILE '{table_file}' INTO TABLE {table_name}
"""
GET_ALL_INDELS = """
SELECT DISTINCT variant_id, POS, REF, ALT
FROM variant_chr{CHROM}
WHERE indel = 1
"""
GET_TRANSLATION_MD5_ID = """
SELECT translation_md5_id
FROM homo_sapiens_variation_74_37.transcript_translation_mapping
WHERE stable_id = "{stable_id}"
"""
GET_POLYPHEN_PREDICTION_MATRIX = """
SELECT prediction_matrix
FROM homo_sapiens_variation_74_37.protein_function_predictions
WHERE translation_md5_id = {translation_md5_id}
    AND analysis_attrib_id = {attrib_id}
"""
GET_SAMPLE_NAME = """
SELECT sample_name
FROM sample
WHERE sample_id = {sample_id}
"""
GET_SAMPLE_ID = """
SELECT sample_id
FROM sample
WHERE sample_name = "{sample_name}" AND sample_type = "{sample_type}"
    AND capture_kit = "{capture_kit}" AND prep_id = {prep_id}
"""
INSERT_SAMPLE = """
INSERT INTO sample (sample_name, sample_type, capture_kit, prep_id)
VALUE ("{sample_name}", "{sample_type}", "{capture_kit}", {prep_id})
"""
CHECK_SAMPLE_IN_SEQDB = """
SELECT 1
FROM seqdbClone
WHERE CHGVID = "{sample_name}" AND SeqType = "{sample_type}" AND
    ExomeSamPrepKit = "{capture_kit}" AND prepId = {prep_id}
"""
