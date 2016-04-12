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
POLYPHEN_ATTRIB_ID = {"humvar":268, "humdiv":269}
# PolyPhen scores are packed by sorted one-letter codes
AMINO_ACIDS = dict(
    [[aa, aa_idx] for aa_idx, aa in enumerate((
        "Ala", #A
        "Cys", #C
        "Asp", #D
        "Glu", #E
        "Phe", #F
        "Gly", #G
        "His", #H
        "Ile", #I
        "Lys", #K
        "Leu", #L
        "Met", #M
        "Asn", #N
        "Pro", #P
        "Gln", #Q
        "Arg", #R
        "Ser", #S
        "Thr", #T
        "Val", #V
        "Trp", #W
        "Tyr" #Y
    ))])
# the regex SnpEff follows for outputting missense changes
HGVS_P_PATTERN = (r"p\.[A-Z][a-z]{2}(?P<codon_position>\d+)"
                  "(?P<amino_acid_change>[A-Z][a-z]{2})")
POLYPHEN_PROB_BITMASK = 2 ** 10 - 1
