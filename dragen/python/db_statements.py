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
