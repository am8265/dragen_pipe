#!/usr/bin/env python
"""
Strip all SnpEff annotations out of a VCF to return it to its "unannotated"
state
"""

import argparse
import sys
from utilities import fh

# indicates SnpEff header for deletion
SNPEFF_HEADERS = (
    "##SnpEff", "##INFO=<ID=ANN,", "##INFO=<ID=LOF,", "##INFO=<ID=NMD,")
SNPEFF_TAGS = ("ANN=", "LOF=", "NMD=")
HEADER = "#CHROM"
INFO_FIELD = "INFO"

def main(vcf_fh, output_fh):
    for line in vcf_fh:
        line = line.rstrip("\r\n")
        skip_line = False
        for snpeff_header in SNPEFF_HEADERS:
            if line.startswith(snpeff_header):
                skip_line = True
                break
        if not skip_line:
            output_fh.write(line + "\n")
            if line.startswith(HEADER):
                break
    info_idx = line.split("\t").index(INFO_FIELD)
    for line in vcf_fh:
        fields = line.rstrip("\r\n").split("\t")
        info_no_snpeff = []
        for info_value in fields[info_idx].split(";"):
            include_annotation = True
            for snpeff_tag in SNPEFF_TAGS:
                if info_value.startswith(snpeff_tag):
                    include_annotation = False
                    break
            if include_annotation:
                info_no_snpeff.append(info_value)
        fields[info_idx] = ";".join(info_no_snpeff)
        output_fh.write("\t".join(fields) + "\n")
    vcf_fh.close()
    if output_fh is not sys.stdout:
        output_fh.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("VCF", type=fh, help="the VCF to process")
    parser.add_argument("-o", "--output", type=argparse.FileType("w"),
                        default=sys.stdout, help="the output VCF")
    args = parser.parse_args()
    main(args.VCF, args.output)
