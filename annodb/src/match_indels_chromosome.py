#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Check whether an indel that's not already in the database exactly is
represented by another one already in it
"""

import MySQLdb
from collections import defaultdict
from pyfaidx import Fasta
from ConfigParser import SafeConfigParser
from db_statements import *
from dragen_globals import *

FLANKING_SIZE = 1000
ALL_INDELS = defaultdict(lambda: defaultdict(list))
indels_queried = False
#nskipped = 0

def match_indel(cur, CHROM, POS, REF, ALT):
    """return the variant_id of a match against the indels currently in the
    database if present, otherwise None
    """
    #global nskipped
    if not indels_queried:
        get_all_indels(cur, CHROM)
    block_id = POS / FLANKING_SIZE
    if (POS - FLANKING_SIZE) < 0 or (POS + FLANKING_SIZE) > chromosome_length:
        # only perform indel matching if the indel is not too close to either
        # end of the chromosome
        return None
    len_REF = len(REF)
    len_indel = len(ALT) - len_REF
    indel_sequence = "".join([
        sequence[POS - FLANKING_SIZE - 1:POS - 1], ALT,
        sequence[POS + len_REF - 1:POS + FLANKING_SIZE - len_indel - 1]])
    for block in xrange(block_id - 1, block_id + 2):
        if block and block in ALL_INDELS and len_indel in ALL_INDELS[block]:
            # only look in the three surrounding blocks where the length of the
            # indel is the same as the candidate novel indel
            for variant_id, db_POS, db_REF, db_ALT in ALL_INDELS[block][len_indel]:
                #if (db_POS == POS and db_REF == REF and db_ALT == ALT):
                #    # debug code to test whether indels match each other in one
                #    # sample
                #    nskipped += 1
                #    continue
                if abs(db_POS - POS) <= FLANKING_SIZE:
                    db_head = sequence[
                        POS - FLANKING_SIZE - 1:db_POS - 1]
                    db_sequence = "".join([
                        db_head, db_ALT,
                        sequence[db_POS + len(db_REF) - 1:db_POS + 2 *
                                 FLANKING_SIZE - len(db_head) - len_indel - 1]])
                    if indel_sequence == db_sequence:
                        return variant_id
    return None

def get_all_indels(cur, CHROM):
    """populate ALL_INDELS with lists of blocks of indels of
    FLANKING_SIZE length for the specified chromosome
    """
    global ALL_INDELS
    config_parser = SafeConfigParser()
    config_parser.read(CNF)
    genome = Fasta(config_parser.get("annodb", "REFERENCE_GENOME"))
    global sequence
    sequence = genome[str(CHROM)][:].seq
    global chromosome_length
    chromosome_length = len(sequence)
    genome.close()
    cur.execute(GET_ALL_INDELS.format(CHROM=CHROM))
    for variant_id, POS, REF, ALT in cur.fetchall():
        len_indel = len(ALT) - len(REF)
        ALL_INDELS[POS / FLANKING_SIZE][len_indel].append(
            (variant_id, POS, REF, ALT))
    global indels_queried
    indels_queried = True
