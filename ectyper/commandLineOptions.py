#!/usr/bin/env python

import argparse
from ectyper import __version__


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
        ivalue = int(value)
        if ivalue <= 0 or ivalue > 100:
            raise argparse.ArgumentTypeError(
                "{0} is an invalid positive int percentage value".format(value)
            )
        return ivalue

    parser = argparse.ArgumentParser(
        description='ectyper v{} Prediction of Escherichia coli serotype from raw reads'
            ' or assembled genome sequences'.format(__version__)
    )

    parser.add_argument(
        "-V",
        "--version",
        action='version',
        version="%(prog)s {}".format(__version__)
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
        help="The number of cores to run ectyper with",
        default=1
    )

    parser.add_argument(
        "-d",
        "--percentIdentity",
        type=check_percentage,
        help="Percent identity required for an allele match [default 90]",
        default=90
    )

    parser.add_argument(
        "-l",
        "--percentLength",
        type=check_percentage,
        help="Percent length required for an allele match [default 50]",
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
        help="Location of pre-computed MASH RefSeq sketch. If provided, genomes "
             "identified as non-E. coli will have their species identified using "
             "MASH. For best results the pre-sketched RefSeq archive "
             "https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh "
             "is recommended"
    )

    if args is None:
        return parser.parse_args()
    else:
        return parser.parse_args(args)
