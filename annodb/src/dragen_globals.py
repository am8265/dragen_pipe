"""
Constant variables shared across modules for the DRAGEN pipeline
"""

CNF = "/nfs/goldstein/software/dragen/dragen.cnf" # defaults file for pipeline
VCF_COLUMNS = ["CHROM", "POS", "rs_number", "REF", "ALT", "QUAL", "FILTER",
               "INFO", "FORMAT", "call"]
BLOCK_SIZE = 10000 # the bases in a block of variant calls (for indexing)
VALID_GTS = set(["0", "1"]) # valid values in GT field, i.e. REF/ALT
# the table format to output for calls
VARIANT_CALL_FORMAT = ("{" + "}\t{".join(
    ["sample_id", "variant_id", "block_id", "GT", "DP_pileup", "DP", "AD_REF",
     "AD_ALT", "GQ", "FS", "MQ", "QD", "QUAL", "ReadPosRankSum", "MQRankSum",
     "PASS"]) + "}")
