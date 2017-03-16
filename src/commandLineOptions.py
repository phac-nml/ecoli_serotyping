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
                        type=int,
                        help="Trigger the use of the E. coli serotyper. \
                              Options are 0 and 1. Default is 0 (false).",
                        default=0,
                        choices=[0, 1])

    parser.add_argument("-f",
                        "--viruelenceFactors",
                        type=int,
                        help="Trigger the use of the Virulence Factors tool. \
                              Options are 0 and 1. Default is 0 (false).",
                        default=0,
                        choices=[0, 1])

    parser.add_argument("-i",
                        "--percentIdentity",
                        type=int,
                        help="Percentage of identity wanted to use against the\
                              database. From 0 to 100, default is 90%.",
                        default=90,
                        choices=range(0, 100))

    parser.add_argument("-l",
                        "--percentLength",
                        type=int,
                        help="Percentage of length wanted to use against the \
                              database. From 0 to 100, default is 90%.",
                        default=90,
                        choices=range(0, 100))

    parser.add_argument("-m",
                        "--minimumGenomes",
                        type=int,
                        help="Minimum number of genomes threshold for a \
                              virulence factor. Default is 1.",
                        default=1)

    parser.add_argument("-c",
                        "--csv",
                        type=int,
                        help="If set to 1, the results will be sent to a .csv \
                              file in the temp/Results folder. Options are 0 \
                              and 1, default=1.",
                        default=1,
                        choices={0, 1})

    return parser.parse_args()
