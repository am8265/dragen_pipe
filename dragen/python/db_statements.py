GET_SAMPLE_INFO = """
SELECT sample_name, sample_type, capture_kit, prep_id
FROM sample
WHERE pseudo_prepid = {pseudo_prepid}
"""
GET_PIPELINE_STEP_ID = """
SELECT id
FROM pipeline_step
WHERE step_name = "{step_name}"
"""
GET_PIPELINE_FINISHED_ID = """
SELECT id
FROM pipeline_step
WHERE step_name = "ArchiveSample""
"""
GET_STEP_STATUS = """
SELECT finished
FROM sample_pipeline_step
WHERE pseudo_prepid = {pseudo_prepid} AND pipeline_step_id = {pipeline_step_id}
"""
INSERT_PIPELINE_STEP = """
INSERT INTO sample_pipeline_step (pseudo_prepid, pipeline_step_id)
VALUE ({pseudo_prepid}, {pipeline_step_id})
"""
UPDATE_PIPELINE_STEP_SUBMIT_TIME = """
UPDATE sample_pipeline_step
SET submit_time = NOW()
WHERE pseudo_prepid = {pseudo_prepid} AND pipeline_step_id = {pipeline_step_id}
"""
UPDATE_PIPELINE_STEP_FINISH_TIME = """
UPDATE sample_pipeline_step
SET finish_time = NOW(), finished = 1
WHERE pseudo_prepid = {pseudo_prepid} AND pipeline_step_id = {pipeline_step_id}
"""
UPDATE_PIPELINE_FINISH_SUBMIT_TIME = """
UPDATE sample_pipeline_step s1
INNER JOIN (
SELECT pseudo_prepid, MIN(submit_time) as st
FROM sample_pipeline_step
WHERE pseudo_prepid = {pseudo_prepid} AND pipeline_step_id <> {pipeline_step_id}
GROUP BY pseudo_prepid) AS s2
ON s1.pseudo_prepid = s2.pseudo_prepid AND s1.pipeline_step_id = {pipeline_step_id}
SET s1.submit_time = s2.st
"""
