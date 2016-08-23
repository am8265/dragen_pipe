GET_SAMPLE_INFO = """
SELECT CHGVID, DBID, prepID, seqtype
FROM seqdbClone
WHERE pseudo_prepid = {pseudo_prepid}
"""
GET_PIPELINE_STEP_ID = """
SELECT id
FROM dragen_pipeline_step_desc
WHERE step_name = "{step_name}"
"""
GET_PIPELINE_FINISHED_ID = """
SELECT id
FROM dragen_pipeline_step_desc
WHERE step_name = "ArchiveSample""
"""
GET_STEP_STATUS = """
SELECT finished
FROM dragen_pipeline_step
WHERE pseudo_prepid = {pseudo_prepid} AND pipeline_step_id = {pipeline_step_id}
"""
INSERT_PIPELINE_STEP = """
INSERT INTO dragen_pipeline_step (pseudo_prepid, pipeline_step_id)
VALUE ({pseudo_prepid}, {pipeline_step_id})
"""
UPDATE_PIPELINE_STEP_SUBMIT_TIME = """
UPDATE dragen_pipeline_step
SET submit_time = NOW()
WHERE pseudo_prepid = {pseudo_prepid} AND pipeline_step_id = {pipeline_step_id}
"""
UPDATE_PIPELINE_STEP_FINISH_TIME = """
UPDATE dragen_pipeline_step
SET finish_time = NOW(), finished = 1
WHERE pseudo_prepid = {pseudo_prepid} AND pipeline_step_id = {pipeline_step_id}
"""
