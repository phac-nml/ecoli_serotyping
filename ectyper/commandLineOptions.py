#!/usr/bin/env python

import argparse

VERSION = "0.3.5"

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
            ' or assembled genome sequences'.format(VERSION)
    )

    parser.add_argument(
        "-V",
        "--version",
        action='version',
        version="%(prog)s {}".format(VERSION)
    )

    parser.add_argument(
        "-i",
        "--input",
        help="Location of E. coli genome file(s). Can be a single file or \
            a directory",
        required=True
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
        "-s",
        "--species",
        action="store_true",
        help="Enable species identification when a non-E. coli genome is found\n\
            Note: refseq downloading is required when running this option for\n\
            the first time"
    )

    parser.add_argument(
        "--detailed",
        action='store_true',
        help="Enable detailed program output"
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Directory location of output files."
    )

    if args is None:
        return parser.parse_args()
    else:
        return parser.parse_args(args)
