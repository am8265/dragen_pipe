#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Take a VCF and annotate with SnpEff
Uses /nfs/goldstein/software/dragen/dragen.cnf to get
optional parameters that are not explicitly overwritten

Brett Copeland <bc2675@cumc.columbia.edu>
3/2016
"""

import argparse
from ConfigParser import SafeConfigParser
from utilities import *

SNPEFF_COMMAND = (
    "{java} -Xmx5G -jar {snpeff} eff {genome_version} -c {snpeff_config} -v "
    "-interval {intervals} -noMotif -noInteraction -noNextProt "
    "-s {output_vcf}.annotations -o vcf {input_vcf} {threaded} > {output_vcf}")
CNF = "/nfs/goldstein/software/dragen/dragen.cnf"

def main(input_vcf, output_vcf, parameters):
    """run the SnpEff command given the specified parameters
    """
    print(SNPEFF_COMMAND.format(
        input_vcf=input_vcf, output_vcf=output_vcf, **parameters))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=CustomFormatter, description=__doc__)
    parser.add_argument("INPUT_VCF", type=file_exists, help="the VCF to annotate")
    parser.add_argument("OUTPUT_VCF", type=create_directory_for_file_if_needed,
                        help="the output VCF file")
    parser.add_argument("--java", type=file_exists,
                        help="custom java program to use")
    parser.add_argument("--snpeff", type=file_exists,
                        help="custom SnpEff jar file to use")
    parser.add_argument("--genome_version", help="custom genome version to use")
    parser.add_argument("--snpeff_config", type=file_exists,
                        help="custom SnpEff configuration to use")
    parser.add_argument("--intervals", type=file_exists,
                        help="custom exome BED file to use")
    parser.add_argument("--single_threaded", default=False, action="store_true",
                        help="only use a single thread")
    args = parser.parse_args()
    config_parser = SafeConfigParser()
    config_parser.read(CNF)
    parameters = {}
    parameters["threaded"] = "" if args.single_threaded else "-t"
    for parameter in (
        "java", "snpeff", "genome_version", "snpeff_config", "intervals"):
        if args.__dict__[parameter]:
            parameters[parameter] = args.__dict__[parameter]
        else:
            parameters[parameter] = config_parser.get("annodb", parameter.upper())
    main(args.INPUT_VCF, args.OUTPUT_VCF, parameters)
