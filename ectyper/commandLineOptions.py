#!/usr/bin/env python

import argparse
from ectyper import __version__
import json, os
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
        description='ectyper v{} antigen database v{}. Prediction of Escherichia coli serotype, pathotype and shiga toxin tying from '
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
        help="Percent identity required for an O antigen allele match [default %(default)s]",
        default=90
    )

    parser.add_argument(
        "-hpid",
        "--percentIdentityHtype",
        type=check_percentage,
        help="Percent identity required for an H antigen allele match [default %(default)s]",
        default=95
    )

    parser.add_argument(
        "-opcov",
        "--percentCoverageOtype",
        type=check_percentage,
        help="Minimum percent coverage required for an O antigen allele match [default %(default)s]",
        default=90
    )

    parser.add_argument(
        "-hpcov",
        "--percentCoverageHtype",
        type=check_percentage,
        help="Minimum percent coverage required for an H antigen allele match [default %(default)s]",
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
        "--reference",
        default=definitions.SPECIES_ID_SKETCH,
        help="Location of pre-computed MASH sketch for species identification. If provided, "
             "genomes "
             "identified as non-E. coli will have their species identified "
             "using "
             "MASH dist"
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
        help="Path to a custom database of O and H antigen alleles in JSON format.\n"
    )

    parser.add_argument(
        "--pathotype",
        action="store_true",
        help="Predict E.coli pathotype and Shiga toxin subtype(s) if present\n"
    )

    parser.add_argument(
        "-pathpid",
        "--percentIdentityPathotype",
        type=check_percentage,
        help="Minimum percent identity required for a pathotype reference allele match [default: %(default)s]",
        default=90
    )

    parser.add_argument(
        "-pathpcov",
        "--percentCoveragePathotype",
        type=check_percentage,
        help="Minimum percent coverage required for a pathotype reference allele match [default: %(default)s]",
        default=50
    )

    if args is None:
        return parser.parse_args()
    else:
        return parser.parse_args(args)
