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
WHERE step_name = "ArchiveSample"
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
SET submit_time = NOW(), times_ran = times_ran + 1
WHERE pseudo_prepid = {pseudo_prepid} AND pipeline_step_id = {pipeline_step_id}
"""
UPDATE_PIPELINE_STEP_FINISH_TIME = """
UPDATE dragen_pipeline_step
SET finish_time = NOW(), finished = 1
WHERE pseudo_prepid = {pseudo_prepid} AND pipeline_step_id = {pipeline_step_id}
"""
GET_DBID_MAX_PREPID = """
SELECT DBID,MAX(p.prepID)
FROM prepT p JOIN pseudo_prepid pp ON p.prepID=pp.prepID WHERE
pseudo_prepID = {pseudo_prepid}
"""
GET_USERID = """
SELECT userid
FROM users WHERE netID = "{userName}"
"""
GET_CAPTURE_KIT_BED = """
SELECT region_file_lsrc
FROM captureKit ck
JOIN prepT p ON ck.prepT_name=p.exomekit
JOIN pseudo_prepid pp ON p.prepid=pp.prepid
WHERE chr="all" and pseudo_prepid={pseudo_prepid}
"""
GET_CAPTURE_KIT = """
SELECT name
FROM captureKit ck
JOIN prepT p ON ck.prepT_name=p.exomekit
JOIN pseudo_prepid pp ON p.prepid=pp.prepid
WHERE chr="all" and pseudo_prepid={pseudo_prepid}
"""
GET_PANEL_OF_NORMALS = """
SELECT file_loc
FROM dragen_panel_of_normals
WHERE capture_kit="{captureKit}" AND sample_type="{sample_type}"
"""
GET_READ_LENGTH = """
"""
GET_SAMPLES_TO_RUN = (
    """ SELECT D.pseudo_prepid, D.sample_name, D.sample_type, D.capture_kit, D.priority """
    """FROM dragen_sample_metadata as D INNER JOIN (select d1.pseudo_prepid FROM """
    """dragen_pipeline_step as d1 WHERE d1.pseudo_prepid NOT IN (SELECT pseudo_prepid """
    """FROM dragen_pipeline_step as d2 WHERE d2.pipeline_step_id = 31 AND d2.finished = 1) """
    """AND d1.pipeline_step_id = 11 AND d1.finished = 1) as P ON """
    """D.pseudo_prepid = P.pseudo_prepid ORDER BY D.priority; """
)
