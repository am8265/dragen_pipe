PROD_GET_SAMPLES = """
SELECT C.CHGVID, C.BioinfoPriority, C.prepid, C.SeqType, C.AlignSeqFileLoc
FROM seqdbClone C
LEFT JOIN ethnicity_status E ON C.CHGVID = E.CHGVID AND C.prepid = E.prepid
    AND C.SeqType = E.SeqType
WHERE C.SeqType IN ("Exome", "Genome") AND (E.predict = 0 OR E.predict IS NULL)
    AND C.AlignSeqFileLoc IS NOT NULL AND C.CHGVID NOT LIKE "CGND%" AND C.Status IN
    ("Completed AnnoDB Pipeline Step Parsing VCF",
    "Failed AnnoDB Pipeline Step Parsing VCF", "In Annotation DB",
    "In Annotation DB/Sample Rejected", "Initialized in AnnoDB - waiting to load",
    "Loaded to Mastery AnnoDB Server - Ready to Replicate", "Passed BioInfo QC",
    "QC done/sample rejected", "QC failure/in Annotation DB", "QC review needed", 
    "Started AnnoDB Pipeline Step VCF Annotation")
"""
UPDATE_PROBS = """
UPDATE seqdbClone SET Caucasian_prob = {}, MiddleEastern_prob = {},
    Hispanic_prob = {}, EastAsian_prob = {}, SouthAsian_prob = {},
    African_prob = {}, genotyping_rate = {}
WHERE prepID = {} AND CHGVID = "{}" AND SeqType = "{}"
"""
UPDATE_PED_STATUS = """
UPDATE ethnicity_status SET create_ped = 1
WHERE prepid = {} AND CHGVID = "{}" AND SeqType = "{}"
"""
RESET_PED_STATUS = """
UPDATE ethnicity_status SET create_ped = 0
WHERE prepid = {} AND CHGVID = "{}" AND SeqType = "{}"
"""
UPDATE_PREDICT_STATUS = """
UPDATE ethnicity_status SET predict = 1
WHERE prepid = {} AND CHGVID = "{}" AND SeqType = "{}"
"""
ETHNICITY_STATUS = """
SELECT {}
FROM ethnicity_status
WHERE prepid = {} AND CHGVID = "{}" AND SeqType = "{}"
"""
INSERT_STATUS = """
INSERT INTO ethnicity_status (prepid, CHGVID, SeqType)
VALUE ({}, "{}", "{}")
"""
GET_FAMILY_ID = """
SELECT FamilyID
FROM SampleT
WHERE CHGVID = "{}"
"""
