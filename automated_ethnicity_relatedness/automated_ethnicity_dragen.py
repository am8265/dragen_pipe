import MySQLdb
from automated_relatedness import get_connection

def update_probs(output_prob_file,pseudo_prepid):
    """
    """
    
    db = get_connection(db="seqdb")
    cmd = (
        """ UPDATE dragen_qc_metrics SET Caucasian_prob = {0},MiddleEastern_prob = {1},"""
        """ Hispanic_prob = {2}, EastAsian_prob = {3} , SouthAsian_prob = {4},"""
        """ African_prob = {5}, genotyping_rate = {6} WHERE pseudo_prepid = {7} """
    )
    
    with open(output_prob_file,'r') as IN:
        for line in IN:
            line = line.strip('\n')
            caucasian,middleeast,hispanic,eastasian,southasian,african,genorate = line.split('\t')
            cur = db.cursor()
            db_cmd = cmd.format(caucasian,middleeast,hispanic,eastasian,southasian,african,genorate,pseudo_prepid)
            print db_cmd
            cur.execute(db_cmd)
            db.commit()
        db.close()
    
        
def get_samples_to_predict():
    """
    Get samples to predict
    """

# WHERE d1.pipeline_step_id = 31 AND d1.finished = 1 ) as O

    query = """
SELECT S.sample_name, S.pseudo_prepid, S.sample_type, S.AlignSeqFileLoc
FROM (SELECT O.pseudo_prepid, O.sample_name, O.sample_type, O.capture_kit, Q.AlignSeqFileLoc, O.priority
      FROM (SELECT D.pseudo_prepid, D.sample_name, D.sample_type, D.capture_kit, D.priority
            FROM dragen_sample_metadata as D
            INNER JOIN dragen_pipeline_step as d1 ON d1.pseudo_prepid = D.pseudo_prepid
            WHERE d1.pipeline_step_id = 31 AND d1.step_status = 'completed' ) as O
      INNER JOIN dragen_qc_metrics as Q ON O.pseudo_prepid = Q.pseudo_prepid) as S
LEFT JOIN dragen_ped_status as P ON S.pseudo_prepid = P.pseudo_prepid
ORDER BY S.priority
"""
    db = get_connection(db="seqdb")
    cur = db.cursor()
    cur.execute(query)
    result = cur.fetchall()
    ##### dragen_ped_status is JUST a pipeline table for each of the stages and separate cols!?!?!?!?!?!?
    ##### presumably populates in one go and then tries to run them in one batch and if anything happens it all breaks.
    ##### may need to wipe --base-output-directory and truncate table!?!?
    max_samples = 100
    # max_samples = 8000
    i = 0 
    for res in result[0:max_samples]:
        i+=1
        yield res
  

