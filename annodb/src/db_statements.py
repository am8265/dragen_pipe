"""
MySQL statements used in the pipeline
"""

# check novelty of variant/get its id
VARIANT_EXISTS_QUERY = """
SELECT variant_id
FROM variant_chr{CHROM}
WHERE POS = {POS} AND REF = "{REF}" AND ALT = "{ALT}"
    AND indel_length = {indel_length}
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
SELECT DISTINCT variant_id, POS, REF, ALT, indel_length
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
GET_DATA_DIRECTORY_FOR_SAMPLE = """
SELECT DragenFileLoc
FROM seqdbClone
WHERE CHGVID = "{sample_name}" AND SeqType = "{sample_type}" AND
    ExomeSamPrepKit = "{capture_kit}" AND prepId = {prep_id}
"""
GET_NUM_CALLS_FOR_SAMPLE = """
SELECT COUNT(DISTINCT(c.variant_id))
FROM called_variant_chr{CHROM} c
INNER JOIN variant_chr{CHROM} v ON c.variant_id = v.variant_id
WHERE c.sample_id = {sample_id}
"""
GET_EFFECT_RANKINGS = """
SELECT *
FROM effect_ranking
"""
GET_SAMPLE_INFO = """
SELECT sample_name, sample_type, capture_kit, prep_id
FROM sample
WHERE sample_id = {sample_id}
"""
GET_PIPELINE_STEP_ID = """
SELECT id
FROM pipeline_step
WHERE description = "Imported Chromosome {chromosome}"
"""
GET_PIPELINE_FINISHED_ID = """
SELECT id
FROM pipeline_step
WHERE description = "Sample Finished"
"""
GET_STEP_STATUS = """
SELECT finished
FROM sample_pipeline_step
WHERE sample_id = {sample_id} AND pipeline_step_id = {pipeline_step_id}
"""
INSERT_PIPELINE_STEP = """
INSERT INTO sample_pipeline_step (sample_id, pipeline_step_id)
VALUE ({sample_id}, {pipeline_step_id})
"""
UPDATE_PIPELINE_STEP_SUBMIT_TIME = """
UPDATE sample_pipeline_step
SET submit_time = NOW()
WHERE sample_id = {sample_id} AND pipeline_step_id = {pipeline_step_id}
"""
UPDATE_PIPELINE_STEP_FINISH_TIME = """
UPDATE sample_pipeline_step
SET finish_time = NOW(), finished = 1
WHERE sample_id = {sample_id} AND pipeline_step_id = {pipeline_step_id}
"""
UPDATE_PIPELINE_FINISH_SUBMIT_TIME = """
UPDATE sample_pipeline_step s1
INNER JOIN (
    SELECT sample_id, MIN(submit_time) as st
    FROM sample_pipeline_step
    WHERE sample_id = {sample_id} AND pipeline_step_id <> {pipeline_step_id}
    GROUP BY sample_id) AS s2
    ON s1.sample_id = s2.sample_id AND s1.pipeline_step_id = {pipeline_step_id}
SET s1.submit_time = s2.st
"""
