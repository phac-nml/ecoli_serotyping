#!/usr/bin/env python

import argparse


def parse_command_line(args=None):
    """
    The options
    for both the serotyper, and virulence finder.
    The returned object is used by both, but the options do not
    necessarily apply to both.
    """

    def check_percentage(value):
        ivalue = int(value)
        if ivalue <= 0 or ivalue > 100:
            raise argparse.ArgumentTypeError("%s is an invalid positive int percentage value" % value)
        return ivalue

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        help="Location of new file(s). Can be a single file or \
            a directory",
        required=True
    )

    parser.add_argument(
        "-d",
        "--percentIdentity",
        type=check_percentage,
        help="Percentage of identity wanted to use against the\
                  database. From 0 to 100, default is 90%%.",
        default=90
    )

    parser.add_argument(
        "-l",
        "--percentLength",
        type=check_percentage,
        help="Percentage of length wanted to use against the \
                  database. From 0 to 100, default is 50%%.",
        default=50
    )

    parser.add_argument(
        "--verify",
        action="store_true",
        help="Enable E. Coli. verification"
    )

    parser.add_argument(
        "-s",
        "--species",
        action="store_true",
        help="Enable non-ecoli species identification"
    )
    
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        help='Enable detailed output'
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Directory location of output files."
    )

    if args is None:
        return parser.parse_args()
    return parser.parse_args(args)
