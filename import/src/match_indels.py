#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Check whether an indel that's not already in the database exactly is
represented by another one already in it
"""

import logging
from collections import defaultdict
from pyfaidx import Fasta
from ConfigParser import SafeConfigParser
from db_statements import *
from waldb_globals import *

def nested_defaultdict(data_type):
    return defaultdict(lambda: defaultdict(data_type))

FLANKING_SIZE = 1000
ALL_INDELS = defaultdict(lambda: nested_defaultdict(list))
chromosome_indels_queried = set()
sequence_by_chromosome = {}
chromosome_lengths = {}
logger = logging.getLogger(__name__)

def get_allele_in_reference_genome(CHROM, POS, REF):
    """get the actual reference allele at the site per the reference genome, for
    either ensuring an added indel is not wrong, or that we don't try to match
    an incorrectly located indel
    """
    return sequence_by_chromosome[CHROM][POS - 1:POS + len(REF) - 1]

def add_new_indel(variant_id, CHROM, POS, REF, ALT, indel_length):
    """add a new indel to the indel set (if it is not incorrectly positioned)
    """
    reference_genome_ref = get_allele_in_reference_genome(CHROM, POS, REF)
    if reference_genome_ref == REF:
        ALL_INDELS[CHROM][POS / FLANKING_SIZE][indel_length].append(
            (variant_id, POS, REF, ALT))
    else:
        logger.warning("REF is wrong for variant {CHROM}-{POS}-{REF}-{ALT}: "
                       "sequence should be {ref}".format(
                           CHROM=CHROM, POS=POS, REF=REF, ALT=ALT,
                           ref=reference_genome_ref))

def match_indel(cur, CHROM, POS, REF, ALT, indel_length):
    """return the variant_id and block_id  of a match against the indels
    currently in the database if present, otherwise None, None
    """
    if CHROM not in chromosome_indels_queried:
        get_all_indels(cur, CHROM)
    # verify the REF allele is correct
    reference_genome_ref = get_allele_in_reference_genome(CHROM, POS, REF)
    if reference_genome_ref != REF:
        logger.warning("REF is wrong for variant {CHROM}-{POS}-{REF}-{ALT}: "
                       "sequence should be {ref}".format(
                           CHROM=CHROM, POS=POS, REF=REF, ALT=ALT,
                           ref=reference_genome_ref))
        return -1, -1

    block_id = POS / FLANKING_SIZE
    if ((POS - FLANKING_SIZE) < 0 or
        (POS + FLANKING_SIZE) > chromosome_lengths[CHROM]):
        # only perform indel matching if the indel is not too close to either
        # end of the chromosome
        return None, None

    len_REF = len(REF)
    indel_sequence = "".join([
        sequence_by_chromosome[CHROM][POS - FLANKING_SIZE - 1:POS - 1], ALT,
        sequence_by_chromosome[CHROM][POS + len_REF - 1:(
            POS + FLANKING_SIZE - indel_length - 1)]])
    for block in xrange(block_id - 1, block_id + 2):
        if (CHROM in ALL_INDELS and block in ALL_INDELS[CHROM] and indel_length in
            ALL_INDELS[CHROM][block]):
            # only look in the three surrounding blocks where the length of the
            # indel is the same as the candidate novel indel
            for variant_id, db_POS, db_REF, db_ALT in \
                 ALL_INDELS[CHROM][block][indel_length]:
                # short-circuit if exact match
                if POS == db_POS and REF == db_REF and ALT == db_ALT:
                    return None, None
                if abs(db_POS - POS) <= FLANKING_SIZE:
                    db_head = sequence_by_chromosome[CHROM][
                        POS - FLANKING_SIZE - 1:db_POS - 1]
                    db_sequence = "".join([
                        db_head, db_ALT,
                        sequence_by_chromosome[CHROM][db_POS + len(db_REF) - 1:(
                            db_POS + 2 * FLANKING_SIZE - len(db_head) -
                            indel_length - 1)]])
                    if indel_sequence == db_sequence:
                        return variant_id, db_POS / BLOCK_SIZE
    return None, None

def get_all_indels(cur, CHROM):
    """populate ALL_INDELS with lists of blocks of indels of
    FLANKING_SIZE length for the specified chromosome
    """
    if CHROM in chromosome_indels_queried:
        logger.warning("{CHROM} was already queried".format(CHROM=CHROM))
        return
    config_parser = SafeConfigParser()
    config_parser.read(CNF)
    genome = Fasta(config_parser.get("annodb", "REFERENCE_GENOME"))
    sequence_by_chromosome[CHROM] = genome[str(CHROM)][:].seq
    chromosome_lengths[CHROM] = len(sequence_by_chromosome[CHROM])
    genome.close()
    cur.execute(GET_ALL_INDELS.format(CHROM=CHROM))
    for variant_id, POS, REF, ALT, indel_length in cur.fetchall():
        add_new_indel(variant_id, CHROM, POS, REF, ALT, indel_length)
    chromosome_indels_queried.add(CHROM)
