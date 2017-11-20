"""
MySQL statements used in the pipeline
"""

# check novelty of variant/get its id
VARIANT_EXISTS_QUERY = """
SELECT variant_id, effect_id, has_high_quality_call
FROM variant_chr{CHROM}
WHERE POS = {POS} AND REF = "{REF}" AND ALT = "{ALT}"
    AND indel_length = {indel_length}
"""
GET_VARIANT_EFFECTS = """
SELECT DISTINCT(effect_id)
FROM variant_chr{CHROM}
WHERE variant_id = {variant_id}
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
FROM homo_sapiens_variation_87_37.transcript_translation_mapping
WHERE stable_id = "{stable_id}"
"""
GET_POLYPHEN_PREDICTION_MATRIX = """
SELECT prediction_matrix
FROM homo_sapiens_variation_87_37.protein_function_predictions
WHERE translation_md5_id = {translation_md5_id}
    AND analysis_attrib_id = {attrib_id}
"""
GET_NUM_CALLS_FOR_SAMPLE = """
SELECT c.variant_id
FROM called_variant_chr{CHROM} c
INNER JOIN variant_chr{CHROM} v ON c.variant_id = v.variant_id
WHERE c.sample_id = {sample_id}
"""
GET_EFFECT_RANKINGS = """
SELECT id, impact, effect
FROM effect_ranking
"""
GET_SAMPLE_INFO = """
SELECT sample_name, sample_type, capture_kit, prep_id
FROM sample
WHERE sample_id = {sample_id}
"""
GET_SAMPLE_INFO_PSEUDO_PREPID = """
SELECT sample_name, sample_type, capture_kit, sample_id
FROM sample
WHERE prep_id = {pseudo_prepid}
"""
GET_DATA_LOADED_PIPELINE_STEP_ID = """
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
SELECT step_status
FROM dragen_pipeline_step
WHERE pseudo_prepid = {prep_id} AND pipeline_step_id = {pipeline_step_id}
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
SET finish_time = NOW(), step_status = "completed"
WHERE pseudo_prepid = {prep_id} AND pipeline_step_id = {pipeline_step_id}
"""
SET_SAMPLE_FINISHED = """
UPDATE sample
SET sample_finished = 1
WHERE sample_id = {sample_id}
"""
INSERT_BIN_STATEMENT = """
LOAD DATA INFILE '{data_file}' IGNORE INTO TABLE {data_type}_bins_chr{chromosome}
    (@dummy, @block_id, @bin_string)
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
GET_SAMPLES_TO_IMPORT = """
SELECT prep_id, sample_name, priority, sample_id, sample_type
FROM sample
WHERE sample_finished = 0{failed_samples_clause}{sample_name_clause}
ORDER BY initialization_time
"""
GET_SAMPLE_DIRECTORY = """
SELECT AlignSeqFileLoc
FROM dragen_qc_metrics
WHERE pseudo_prepid = {prep_id}
"""
STEP_FINISHED = """
SELECT 1
FROM dragen_pipeline_step p
INNER JOIN dragen_pipeline_step_desc d ON p.pipeline_step_id = d.id
WHERE p.pseudo_prepid = {prep_id} AND d.step_name = "Imported {step_name}"
    AND p.step_status = "completed"
"""
GET_PIPELINE_STEP_ID = """
SELECT id
FROM dragen_pipeline_step_desc
WHERE step_name = "{step_name}"
"""
GET_SAMPLES_TO_INITIALIZE = """
SELECT p1.pseudo_prepid
FROM dragen_pipeline_step p1
LEFT JOIN dragen_pipeline_step p2 ON p1.pseudo_prepid = p2.pseudo_prepid
    AND p1.pipeline_step_id = {sample_archived_step_id} AND p1.step_status = "completed"
    AND p2.pipeline_step_id = {sample_initialized_step_id} AND p2.step_status = "completed"
WHERE p1.pipeline_step_id = {sample_archived_step_id} AND p1.step_status = "completed"
    AND p2.pseudo_prepid IS NULL
"""
GET_SAMPLE_METADATA = """
SELECT sample_name, sample_type, capture_kit, seqscratch_drive, priority
FROM dragen_sample_metadata
WHERE pseudo_prepid = {prep_id}
"""
GET_CAPTURE_KIT_BED = """
SELECT region_file_lsrc
FROM captureKit
WHERE chr = "all" AND name = "{capture_kit}"
"""
INITIALIZE_SAMPLE = """
INSERT INTO sample (sample_name, sample_type, capture_kit, prep_id, priority)
VALUE ("{sample_name}", "{sample_type}", "{capture_kit}", {prep_id}, {priority})
"""
INITIALIZE_SAMPLE_SEQDB = """
REPLACE dragen_pipeline_step (pseudo_prepid, pipeline_step_id,
    version, finish_time, times_ran, step_status)
VALUE ({pseudo_prepid}, {sample_initialized_step_id}, "{version}", NOW(), 1, "completed")
"""
BEGIN_STEP = """
REPLACE dragen_pipeline_step
    (pseudo_prepid, pipeline_step_id, version, submit_time, finish_time, times_ran, step_status)
VALUE ({pseudo_prepid}, {pipeline_step_id}, "{version}", NOW(), NULL, {times_ran}, "started")
"""
FINISH_STEP = """
UPDATE dragen_pipeline_step
SET step_status = "completed", finish_time = NOW()
WHERE pseudo_prepid = {pseudo_prepid} AND pipeline_step_id = {pipeline_step_id}
"""
FAIL_STEP = """
UPDATE dragen_pipeline_step
SET step_status = "failed", finish_time = NOW()
WHERE pseudo_prepid = {pseudo_prepid} AND pipeline_step_id = {pipeline_step_id}
"""
