#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Check the integrity of a BAM from the pipeline and return any errors
Checks:
    proper compression
    minimum number of reads per chromosome (optional, may not be appropriate for custom capture)
    proper read groups
"""

import pysam
import sys
import subprocess
import os
from collections import Counter
CHROMOSOMES = [str(c) for c in range(1, 23)] + ["X"]
CHROMOSOMES_SET = set(CHROMOSOMES)
# should obviously be a stricter test than this, but for now just check each
# chromosome has a read to detect truncation issues
MINIMUM_READS_PER_CHROMOSOME = 1

def check_bam(bam_fn, check_read_counts=True):
    if not os.path.isfile(bam_fn):
        raise OSError("{} does not exist!".format(bam_fn))
    bam = pysam.AlignmentFile(bam_fn, "rb")
    errors = []
    read_group_mismatch_detected = False
    bam_corruption_detected = False
    read_group_missing_detected = False
    extra_read_group_detected = False
    try:
        if check_read_counts:
            read_counts_by_chromosome = Counter()
        read_groups = {}
        for read_group in bam.header["RG"]:
            # the ID in our BAMs is the flowcell.lane, but the FASTQ record will
            # be flowcell:lane
            read_groups[read_group["ID"]] = ":" + read_group["PU"].replace(".", ":") + ":"
        for x, read in enumerate(bam.fetch(until_eof=True)):
            if check_read_counts:
                if not read.is_unmapped:
                    chromosome = read.reference_name
                    if chromosome in CHROMOSOMES_SET:
                        read_counts_by_chromosome[chromosome] += 1
            if read.has_tag("RG"):
                read_group = read.get_tag("RG")
            else:
                if not read_group_missing_detected:
                    errors.append("Read #{} is missing a read group tag.".format(
                        x + 1))
                    read_group_missing_detected = True
                continue
            if read_group not in read_groups:
                if not extra_read_group_detected:
                    errors.append("Read #{} refers to a read group ({}) that is not "
                                  "present in the header.".format(x + 1, read_group))
                    extra_read_group_detected = True
                continue
            if read_groups[read_group] not in read.query_name:
                if not read_group_mismatch_detected:
                    errors.append("Read #{}'s QNAME ({}) does not match the PU + "
                                  "lane in the header for RG {}: {}".format(
                                      x + 1, read.query_name, read_group,
                                      read_groups[read_group]))
                    read_group_mismatch_detected = True
                continue
        if check_read_counts:
            for chromosome in CHROMOSOMES:
                if read_counts_by_chromosome[chromosome] < MINIMUM_READS_PER_CHROMOSOME:
                    errors.append("Chromosome {} does not meet the minimal "
                                  "threshold of reads ({}): {}!".format(
                                      chromosome, MINIMUM_READS_PER_CHROMOSOME,
                                      read_counts_by_chromosome[chromosome]))
    finally:
        if not bam.closed:
            bam.close()
    with open(os.devnull, "w") as devnull:
        p = subprocess.Popen(["samtools", "view", bam_fn], stdout=devnull,
                             stderr=subprocess.PIPE)
    _, err = p.communicate()
    if p.returncode or "EOF marker is absent" in err:
        errors.append("BAM file is corrupted!")
    return errors


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("BAM", help="the BAM file to check")
    parser.add_argument("-c", "--check_read_counts", default=False,
                        action="store_true", help="confirm whether the BAM has a "
                        "minimal threshold of reads across each chromosome "
                        "(generally should be used except for custom capture " 
                        "samples)")
    args = parser.parse_args()
    errors = check_bam(args.BAM, args.check_read_counts)
    for error in errors:
        sys.stderr.write(error + "\n")
        sys.exit(1)
