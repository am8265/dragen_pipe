import MySQLdb
from automated_relatedness import get_connection

def update_probs(output_prob_file,prepid,sample_name,sample_type):
    """
    """
    
    db = get_connection(db="seqdb")
    cmd = (
        """ UPDATE seqdbClone SET Caucasian_prob = {0},MiddleEastern_prob = {1},"""
        """ Hispanic_prob = {2}, EastAsian_prob = {3} , SouthAsian_prob = {4},"""
        """ African_prob = {5}, genotyping_rate = {6} WHERE prepID = {7} AND"""
        """ CHGVID = '{8}' AND SeqType = '{9}'"""
    )
    
    with open(output_prob_file,'r') as IN:
        for line in IN:
            line = line.strip('\n')
            caucasian,middleeast,hispanic,eastasian,southasian,african,genorate = line.split('\t')
            cur = db.cursor()
            db_cmd = cmd.format(caucasian,middleeast,hispanic,eastasian,southasian,african,genorate,
                                   prepid,sample_name,sample_type)
            print db_cmd
            cur.execute(db_cmd)
            db.commit()
        db.close()
    
        
def get_samples_to_predict():
    """
    Get samples to predict
    """

    query = ("""
             SELECT S.CHGVID, S.prepid, S.SeqType, S.AlignSeqFileLoc, S.BioinfoPriority
             FROM seqdbClone as S
             LEFT JOIN ethnicity_status as E ON S.CHGVID = E.CHGVID
                 AND S.prepid = E.prepid AND S.SeqType = E.SeqType
             WHERE S.SeqType IN ('Exome','Genome') AND S.Status LIKE '%Annotation DB%'
                AND (E.predict = 0 OR E.predict IS NULL) AND S.AlignSeqFileLoc IS NOT NULL
             ORDER BY S.BioinfoPriority"""
        )
    db = get_connection(db="seqdb")
    cur = db.cursor()
    cur.execute(query)
    result = cur.fetchall()
    max_samples = 10000
    i = 0 
    for res in result[0:max_samples]:
        i+=1
        yield res
  
def update_ped_status(chgvid,prepid,seqtype,field):
    """ Update status
    chgvid : str
    prepid : str
    seqtype : str
    field : str 
    """

    update_statement = (""" UPDATE ethnicity_status SET {0} = 1"""
                        """ WHERE prepid = {1} AND seqType = '{2}' AND"""
                        """ CHGVID = '{3}' """.format(
                            field,prepid,seqtype,chgvid)
                        )            
    db = get_connection(db="seqdb")
    try:
        cur = db.cursor()
        cur.execute(update_statement)
        db.commit()
    finally:
        db.close()


