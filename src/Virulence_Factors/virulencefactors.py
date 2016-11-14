#!/usr/bin/env python

import logging
import argparse
import os
import subprocess

from src.createdirs import createDirs

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"

def parseCommandLine():
    """
    Initalizing the two main commands of the command line for the project.
    - input: refers to the location of the file(s) that will be processed
    - out: refers to the output of the program. Default is STDOUT
    - pi: refers to the percentage of identity wanted. Default is 90%
    - pl: refers to the percentage of length wanted. Default is 90%.
    - v: refers to the verbosity.
    - csv: if the user wants a csv copy of the results.

    :return parser.parse_args(): Data from the commands.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-in", "--input", help="Location of new file(s). Can be a single file or a directory")
    parser.add_argument("-out", "--output", type=argparse.FileType('w'), help="Output of the program. Default is STDOUT.", default=sys.stdout)
    parser.add_argument("-pi", "--percentIdentity", type=int, help="Percentage of identity wanted to use against the database. From 0 to 100, default is 90%.", default=90)
    parser.add_argument("-pl", "--percentLength", type=int, help="Percentage of length wanted to use against the database. From 0 to 100, default is 90%.", default=90)
    parser.add_argument("-csv", help="If set to true, the results will be sent to a .csv file in the temp/Results folder.", default='true')

    return parser.parse_args()


def initializeDB():
    """
    Generating the database if it was not already generated. The database can be found in the tmp folder.

    :return int 0 or 1: 0 being that the database was created successfully (or already existed).
    """

    REL_DIR = SCRIPT_DIRECTORY + '../temp/databases/VF_Database/'

    if not os.path.isdir(REL_DIR):
        os.mkdir(REL_DIR)

    if os.path.isfile(REL_DIR + 'VirulenceFactorsDB.nin'):
        return 0
    else:
        return subprocess.call(["/usr/bin/makeblastdb", "-in", SCRIPT_DIRECTORY + "../Data/repaired_ecoli_vfs.ffn ", "-dbtype", "nucl", "-title", "VirulenceFactorsDB", "-out", REL_DIR + "VirulenceFactorsDB"])




if __name__=='__main__':

    logging.basicConfig(filename='ectyper.log',level=logging.INFO)
    args = parseCommandLine()
    createDirs()
