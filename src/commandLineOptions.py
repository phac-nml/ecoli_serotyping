#!/usr/bin/env python

import argparse


def parse_command_line():
    """
        The options for both the serotyper, and virulence finder.
        The returned object is used by both, but the options do not
        necessarily apply to both.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-i",
                        "--input",
                        help="Location of new file(s). Can be a single file or \
                              a directory",
                        required=True)

    parser.add_argument("-s",
                        "--serotyper",
                        action="store_true",
                        help="Trigger the use of the E. coli serotyper.")

    parser.add_argument("-f",
                        "--virulenceFactors",
                        action="store_true",
                        help="Trigger the use of the Virulence Factors tool.")

    parser.add_argument("-d",
                        "--percentIdentity",
                        type=int,
                        help="Percentage of identity wanted to use against the\
                              database. From 0 to 100, default is 90%.",
                        default=90,
                        choices=range(1, 101))

    parser.add_argument("-l",
                        "--percentLength",
                        type=int,
                        help="Percentage of length wanted to use against the \
                              database. From 0 to 100, default is 90%.",
                        default=90,
                        choices=range(1, 101))

    parser.add_argument("-m",
                        "--minimumGenomes",
                        type=int,
                        help="Minimum number of genomes threshold for a \
                              virulence factor. Default is 1.",
                        default=1)

    parser.add_argument("-t",
                        "--tabular",
                        action="store_true",
                        help="Sets the output format to `tabular`. JSON otherwise.")

    return parser.parse_args()
