GET_ARCHIVED_SAMPLES=(
    """ SELECT D.pseudo_prepid,S.CHGVID,S.SeqType,S.capture_kit FROM
    dragen_pipeline_step as D INNER JOIN dragen_sample_metadata as S
    ON D.pseudo_prepid = S.pseudo_prepid AND D.pipeline_step_id = 31
    AND D.finished = 1
    """
)

## Samples to append to master ped
GET_SAMPLES_FOR_APPEND=(
    """ SELECT Q.

