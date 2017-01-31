"""
MySQL statements used in the pipeline
"""

# check novelty of variant/get its id
VARIANT_EXISTS_QUERY = """
SELECT DISTINCT variant_id, effect_id, has_high_quality_call
FROM variant_chr{CHROM}
WHERE POS = {POS} AND REF = "{REF}" AND ALT = "{ALT}"
    AND indel_length = {indel_length}
"""
GET_VARIANT_EFFECTS = """
SELECT DISTINCT(effect_id)
FROM variant_chr{CHROM}
WHERE variant_id = {variant_id}
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
LOAD_TABLE_REPLACE = """
LOAD DATA INFILE '{table_file}' REPLACE INTO TABLE {table_name}
"""
GET_ALL_INDELS = """
SELECT variant_id, POS, REF, ALT, indel_length
FROM indel_chr{CHROM}
"""
MATCHED_INDEL_EXISTS = """
SELECT variant_id
FROM matched_indels
WHERE CHROM = "{CHROM}" AND POS = {POS} AND REF = "{REF}" AND ALT = "{ALT}"
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
    AND capture_kit = "{capture_kit}"
"""
INSERT_SAMPLE = """
INSERT INTO sample (sample_name, sample_type, capture_kit, prep_id, priority)
VALUE ("{sample_name}", "{sample_type}", "{capture_kit}", {prep_id}, {priority})
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
GET_PIPELINE_SAMPLE_INITIALIZED_ID = """
SELECT id
FROM dragen_pipeline_step_desc
WHERE step_name = "Sample Initialized in DB"
"""
GET_PIPELINE_STEP_ID = """
SELECT id
FROM dragen_pipeline_step_desc
WHERE step_name = "Imported Chromosome {chromosome} {data_type}"
"""
GET_PIPELINE_FINISHED_ID = """
SELECT id
FROM dragen_pipeline_step_desc
WHERE step_name = "Sample Finished"
"""
GET_STEP_STATUS = """
SELECT finished
FROM dragen_pipeline_step
WHERE pseudo_prepid = {prep_id} AND pipeline_step_id = {pipeline_step_id}
"""
INITIALIZE_SAMPLE_PIPELINE_STEP = """
INSERT INTO dragen_pipeline_step 
    (pseudo_prepid, pipeline_step_id, finish_time, times_ran, finished)
VALUE ({prep_id}, {pipeline_step_id}, NOW(), 1, 1)
"""
INSERT_PIPELINE_STEP = """
INSERT INTO dragen_pipeline_step (pseudo_prepid, pipeline_step_id)
VALUE ({prep_id}, {pipeline_step_id})
"""
GET_TIMES_STEP_RUN = """
SELECT times_ran
FROM dragen_pipeline_step
WHERE pseudo_prepid = {prep_id} AND pipeline_step_id = {pipeline_step_id}
"""
INCREMENT_TIMES_STEP_RUN = """
UPDATE dragen_pipeline_step
SET times_ran = times_ran + 1
WHERE pseudo_prepid = {prep_id} AND pipeline_step_id = {pipeline_step_id}
"""
UPDATE_PIPELINE_STEP_SUBMIT_TIME = """
UPDATE dragen_pipeline_step
SET submit_time = NOW(), times_ran = {times_run}
WHERE pseudo_prepid = {prep_id} AND pipeline_step_id = {pipeline_step_id}
"""
UPDATE_PIPELINE_STEP_FINISH_TIME = """
UPDATE dragen_pipeline_step
SET finish_time = NOW(), finished = 1
WHERE pseudo_prepid = {prep_id} AND pipeline_step_id = {pipeline_step_id}
"""
UPDATE_PIPELINE_FINISH_SUBMIT_TIME = """
UPDATE dragen_pipeline_step p1
INNER JOIN (
    SELECT pseudo_prepid, MIN(submit_time) as st
    FROM dragen_pipeline_step
    WHERE pseudo_prepid = {prep_id} AND pipeline_step_id NOT IN
        ({pipeline_step_id}, {sample_initialized_id})
    GROUP BY pseudo_prepid) AS p2
ON p1.pseudo_prepid = p2.pseudo_prepid AND p1.pipeline_step_id = {pipeline_step_id}
SET p1.submit_time = p2.st
"""
SET_SAMPLE_FINISHED = """
UPDATE sample
SET sample_finished = 1
WHERE sample_id = {sample_id}
"""
INSERT_BIN_STATEMENT = """
LOAD DATA INFILE '{data_file}' INTO TABLE {data_type}_bins_chr{chromosome}
    (@block_id, @bin_string)
    SET sample_id={sample_id}, block_id=@block_id, {data_type}_string=@bin_string
"""
GET_MIN_CUSTOM_TRANSCRIPT_ID = """
SELECT MIN(id)
FROM custom_transcript_ids_chr{CHROM}
"""
GET_CUSTOM_TRANSCRIPT_IDS = """
SELECT id, transcript_ids
FROM custom_transcript_ids_chr{CHROM}
"""
GET_SAMPLE_PREPID = """
SELECT pseudo_prepid, prepid, BioInfoPriority
FROM seqdbClone
WHERE CHGVID = "{sample_name}" AND SeqType = "{sample_type}"
    AND ExomeSamPrepKit = "{capture_kit}" AND Status <> "Blacklisted"
"""
GET_PREPID = """
SELECT pseudo_prepid
FROM pseudo_prepid
WHERE prepid = {prepid}
"""
GET_SAMPLES_TO_IMPORT = """
SELECT sample_name, priority, sample_id, sample_type
FROM sample
WHERE sample_finished = 0{failed_samples_clause}
ORDER BY initialization_time
"""
GET_SAMPLE_INITIALIZED_STEP_ID = """
SELECT id
FROM dragen_pipeline_step_desc
WHERE step_name = "Sample Initialized In DB"
"""
GET_SAMPLE_DIRECTORY = """
SELECT AlignSeqFileLoc
FROM dragen_qc_metrics
WHERE pseudo_prepid = {prep_id}
"""
