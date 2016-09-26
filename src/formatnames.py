#!/usr/bin/env python

import argparse
import os
import re
from ecvalidatingfiles import *

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"

def parseCommandLine():
    """

    :return parser.parse_args():
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Location of new file(s). Can be a single file or a directory")

    return parser.parse_args()


def formatNames(filesList):

    filesDict = {}

    for file in filesList:
        filesDict[file] = [line.rstrip('\n') for line in open(file)]


    for file, lines in filesDict.iteritems():
        filename = os.path.basename(file)
        with open(SCRIPT_DIRECTORY+ "../temp/NEW_"+filename, 'a+') as f:
            newLines = []
            for line in lines:
                matches = re.findall('(O\d+\w*|H\d+\w*)', line)
                if matches != None:
                    for match in matches:
                        line  = re.sub(r'(O\d+\w*|H\d+\w*|[/]*)', '', line)
                        match = str(match)
                        f.write(line+match+"\n")
                        print f
                        newLines.append(line + match)
            filesDict[file] = newLines





if __name__=='__main__':

    args = parseCommandLine()
    filesList = getFilesList(args.input)
    formatNames(filesList)

