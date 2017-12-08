GET_SAMPLES = """
SELECT m.pseudo_prepid, m.sample_name, m.priority, m.sample_type, m.capture_kit
FROM dragen_sample_metadata m
INNER JOIN dragen_pipeline_step p1 ON m.pseudo_prepid = p1.pseudo_prepid
WHERE p1.pipeline_step_id = 1 AND p1.step_status = "completed"
    AND p1.pseudo_prepid < 8388600{sample_type_clause}
    AND NOT EXISTS (
    SELECT 1
    FROM dragen_pipeline_step p2
    WHERE p1.pseudo_prepid = p2.pseudo_prepid AND p2.pipeline_step_id = 31
        AND p2.step_status = "completed")
"""
GET_STEP_NAMES = """
SELECT d2.step_name
FROM dragen_pipeline_step_desc d1
INNER JOIN dragen_pipeline_step_desc d2 ON d1.id < d2.id
INNER JOIN dragen_pipeline_step_desc d3 ON d2.id <= d3.id
WHERE d1.step_name = "DragenAlignment" AND d3.step_name = "ArchiveSample"
"""
GET_DP_BLOCKS_FILE = """
SELECT path
FROM blocks_by_kit
WHERE capture_kit = "{capture_kit}"
"""
ANY_STEP_FAILED = """
SELECT 1
FROM dragen_pipeline_step
WHERE step_status = "failed" AND pipeline_step_id >= 2 AND pipeline_step_id <= 31
    AND pseudo_prepid = {pseudo_prepid}
"""
