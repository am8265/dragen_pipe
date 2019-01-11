#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Check the integrity of a VCF from the pipeline and return any errors
Currently checks only for non-printable characters and optionally meeting a
minimum threshold of variant counts on the non-contig chromosomes + mitochondria
"""

import sys
import os
from string import printable
from collections import defaultdict

CCDS_BASES = "/nfs/seqscratch_ssd/PIPELINE_DATA/ccds_bases.p"
CHROMOSOMES = [str(c) for c in range(1, 23)] + ["X"]

""" These counts represent intervals of the observed means +- 8 standard
deviations of counts from 1000 samples chosen at random.
None of the samples actually failed any of these criteria, the primary
purpose here is to flag any obviously wrong samples (e.g. file truncation/high
levels of contamination, etc.).
"""
MIN_CCDS_VARIANTS = 9199
MAX_CCDS_VARIANTS = 34561
CCDS_VARIANTS_RANGE_BY_CHROMOSOME = {
    "1":(927, 3597), "2":(542, 2318), "3":(402, 1829), "4":(362, 1344),
    "5":(330, 1617), "6":(330, 2552), "7":(156, 1856), "8":(189, 1177),
    # "9":(267, 1407), "10":(389, 1434), "11":(719, 2303), "12":(510, 1643),
    "9":(267, 1407), "10":(389, 1434), "11":(719, 2303), "12":(510, 1695),
    "13":(152, 634), "14":(181, 1252), "15":(213, 1181), "16":(232, 1550),
    # "17":(451, 2053), "18":(138, 616), "19":(381, 3236), "20":(155, 929),
    "17":(451, 2100), "18":(138, 616), "19":(381, 3236), "20":(155, 929),
    "21":(59, 528), "22":(87, 872), "X":(59, 1028)}

def check_vcf(vcf_fn, check_variant_counts=False, debug=False):
    if not os.path.isfile(vcf_fn):
        raise OSError("{} does not exist!".format(vcf_fn))
    if vcf_fn.endswith(".gz"):
        from gzip import open as gopen
        open_function = gopen
    else:
        open_function = open
    errors = []
    printable_set = set(printable)
    with open_function(vcf_fn) as vcf_fh:
        if check_variant_counts:
            ccds_variant_counts_by_chromosome = defaultdict(int)
            from cPickle import load
            ccds_bases = load(open(CCDS_BASES, "rb"))
            ccds_variant_count = 0
        header_finished = False
        try:
            for x, line in enumerate(vcf_fh):
                if set(line) - printable_set:
                    errors.append("Non-printable character(s) on line {}".format(
                        x + 1))
                    continue
                if not header_finished:
                    if line.startswith("#CHROM"):
                        header_finished = True
                    continue
                if check_variant_counts:
                    fields = line.split("\t")
                    chromosome = fields[0]
                    if chromosome in CHROMOSOMES:
                        position = int(fields[1])
                        if position in ccds_bases[chromosome]:
                            ccds_variant_count += 1
                            ccds_variant_counts_by_chromosome[chromosome] += 1
        except IOError:
            errors.append("Reached line #{} after which there is a gzip error".format(
                x + 1))
        finally:
            if not vcf_fh.closed:
                vcf_fh.close()

        if check_variant_counts:
            if ccds_variant_count < MIN_CCDS_VARIANTS:
                msg = ("Minimum # of CCDS variants ({}) not met: {}".
                       format(MIN_CCDS_VARIANTS, ccds_variant_count))
                errors.append(msg)
                if debug:
                    sys.stderr.write(msg + "\n")
            elif ccds_variant_count > MAX_CCDS_VARIANTS:
                msg = ("Maximum # of CCDS variants ({}) exceeded: {}".
                       format(MAX_CCDS_VARIANTS, ccds_variant_count))
                errors.append(msg)
                if debug:
                    sys.stderr.write(msg + "\n")
            for chromosome in CHROMOSOMES:
                min_value, max_value = CCDS_VARIANTS_RANGE_BY_CHROMOSOME[chromosome]
                msg = ("Variant count on chromosome {} {{}}in permissible "
                       "range [{}, {}]: {}".format(
                           chromosome, min_value, max_value,
                           ccds_variant_counts_by_chromosome[chromosome]))
                if (min_value <= ccds_variant_counts_by_chromosome[chromosome] <=
                    max_value):
                    if debug:
                        sys.stderr.write(msg.format("") + "\n")
                else:
                    msg = msg.format("not ")
                    errors.append(msg)
                    if debug:
                        sys.stderr.write(msg + "\n")
    if debug:
        if errors:
            sys.stderr.write("Tests FAILED!\n")
        else:
            sys.stderr.write("Tests passed.\n")

    return errors

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("VCF", help="the VCF file name")
    parser.add_argument("-c", "--check_variant_counts", default=False,
                        action="store_true", help="confirm whether the VCF's "
                        "variant counts by chromosome and in the CCDS meet "
                        "the minimal standards (generally should be used "
                        "except for custom capture samples)")
    parser.add_argument("-d", "--debug", default=False, action="store_true",
                        help="debug mode (print some extra messages)")
    args = parser.parse_args()
    errors = check_vcf(args.VCF, args.check_variant_counts, args.debug)
    if errors:
        for message in errors:
            sys.stderr.write(message + "\n")
        sys.exit(1)
