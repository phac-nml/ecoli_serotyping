#!/usr/bin/env python

import argparse
import os
import re
from ecvalidatingfiles import *
from Bio import SeqIO

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
TITLES = {}


def parseCommandLine():
    """

    :return parser.parse_args():
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("names", help="")
    parser.add_argument("sequences", help="")

    return parser.parse_args()


def formatNames(filesList):

    filesDict = {}

    for file in filesList:
        filesDict[file] = [line.rstrip('\n') for line in open(file)]


    for file, lines in filesDict.iteritems():
        for line in lines:
            matches = re.findall('(O\d+\w*|H\d+\w*)', line)
            if matches != None:
                oldLine = line
                TITLES[oldLine] = []
                for match in matches:
                    line  = re.sub(r'(O\d+\w*|H\d+\w*|[/]*)', '', line)
                    match = str(match)
                    TITLES[oldLine].append(line+match)




def addSeqToDB(filesList):

    with open(SCRIPT_DIRECTORY + "../Data/EcOH.fasta", 'a+') as handle:
        for file in filesList:
            for record in SeqIO.parse(file, "fasta"):

                if record.id in TITLES:
                    for title in TITLES[record.id]:
                        record.id = title
                        record.description = title
                        record.name = title
                        SeqIO.write(record, handle, "fasta")





if __name__=='__main__':

    args = parseCommandLine()
    namesList = getFilesList(args.names)
    seqsList = getFilesList(args.sequences)
    formatNames(namesList)
    addSeqToDB(seqsList)




