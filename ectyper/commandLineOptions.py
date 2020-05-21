#!/usr/bin/env python

import argparse
from ectyper import __version__
import json
import ectyper.definitions as definitions

def parse_command_line(args=None):
    """
    Options for E. coli serotype prediction.

    Args:
        args: Optional args to be passed to argparse.parse_args()

    Returns:
        The populated argparse Namespace
    """

    def check_percentage(value):
        """
        type checker for percentage input
        """
        try:
            ivalue = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError(
                "{0} is not an integer. A valid positive integer percentage value is required".format(value)
            )
        if ivalue <= 0 or ivalue > 100:
            raise argparse.ArgumentTypeError(
                "{0} is an invalid positive integer percentage value".format(value)
            )
        return ivalue

    def checkdbversion():
        with open(file=definitions.SEROTYPE_ALLELE_JSON) as fp:
            database = json.load(fp)
        return database["version"]

    dbversion = checkdbversion()

    parser = argparse.ArgumentParser(
        description='ectyper v{} database v{} Prediction of Escherichia coli serotype from '
                    'raw reads'
                    ' or assembled genome sequences. The default settings are recommended.'.format(__version__, dbversion)
    )

    parser.add_argument(
        "-V",
        "--version",
        action='version',
        version="%(prog)s {} running database version {}".format(__version__,dbversion)
    )

    parser.add_argument(
        "-i",
        "--input",
        help="Location of E. coli genome file(s). Can be a single file, a \
            comma-separated list of files, or a directory",
        required=True
    )

    parser.add_argument(
        "-c",
        "--cores",
        type=int,
        help="The number of cores to run ectyper with",
        default=1
    )

    parser.add_argument(
        "-opid",
        "--percentIdentityOtype",
        type=check_percentage,
        help="Percent identity required for an O antigen allele match [default 90]",
        default=90
    )

    parser.add_argument(
        "-hpid",
        "--percentIdentityHtype",
        type=check_percentage,
        help="Percent identity required for an H antigen allele match [default 95]",
        default=95
    )

    parser.add_argument(
        "-opcov",
        "--percentCoverageOtype",
        type=check_percentage,
        help="Minumum percent coverage required for an O antigen allele match [default 95]",
        default=90
    )

    parser.add_argument(
        "-hpcov",
        "--percentCoverageHtype",
        type=check_percentage,
        help="Minumum percent coverage required for an H antigen allele match [default 50]",
        default=50
    )

    parser.add_argument(
        "--verify",
        action="store_true",
        help="Enable E. coli species verification"
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Directory location of output files"
    )

    parser.add_argument(
        "-r",
        "--refseq",
        help="Location of pre-computed MASH RefSeq sketch. If provided, "
             "genomes "
             "identified as non-E. coli will have their species identified "
             "using "
             "MASH. For best results the pre-sketched RefSeq archive "
             "https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh "
             "is recommended"
    )

    parser.add_argument(
        "-s",
        "--sequence",
        action="store_true",
        help="Prints the allele sequences if enabled as the final columns of "
             "the output"
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print more detailed log including debug messages"
    )

    parser.add_argument(
        "--dbpath",
        help="Path to a custom database of O and H antigen alleles in JSON format.\nCheck Data/ectyper_database.json for more information"
    )

    if args is None:
        return parser.parse_args()
    else:
        return parser.parse_args(args)
