#!/usr/bin/env python

import json
import sys
import os

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
sys.path.append(os.path.abspath(SCRIPT_DIRECTORY + '../Serotyper/'))
sys.path.append(os.path.abspath(SCRIPT_DIRECTORY + '../'))

from ecvalidatingfiles import *
from createdirs import *


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
    parser.add_argument("-pi", "--percentIdentity", type=int, help="Percentage of idencdtity wanted to use against the database. From 0 to 100, default is 90%.", default=90)
    parser.add_argument("-pl", "--percentLength", type=int, help="Percentage of length wanted to use against the database. From 0 to 100, default is 90%.", default=90)
    parser.add_argument("-csv", help="If set to 1, the results will be sent to a .csv file in the temp/Results folder. Options are 0 and 1, default=1.", default=1)

    return parser.parse_args()

def getResults(genomesList, RGIpath):
    resultDict = {}

    for genome_name in genomesList:
        filename = os.path.basename(genome_name)
        filename = os.path.splitext(filename)

        out = 'temp/Results/' + filename[0]
        print out
        temp_result = subprocess.call(['python', RGIpath + 'rgi.py', "-i", genome_name, "-o", out])

        formatted_out ='temp/Results/' + filename[0] + 'FORMATTED'
        temp_result = subprocess.call(['python',  RGIpath + 'formatJson.py', "-i", SCRIPT_DIRECTORY + '../../' + out + '.json', "-o", formatted_out])

        with open(SCRIPT_DIRECTORY + '../../' + formatted_out + '.json', 'r') as genome_file:
            resultDict[genome_name] = json.load(genome_file)

    return resultDict



if __name__ == '__main__':

    logging.basicConfig(filename=SCRIPT_DIRECTORY + 'rgi.log',level=logging.INFO)
    logging.info("This log file doesn't contain the logs from the actual RGI tool. To see these logs, go to your downloaded RGI folder and open the app.log file.")

    RGIpath = ''
    for relpath in sys.path:
        if 'release-rgi-v3.1.1-58cad6a3b443abb290cf3df438fe558bc5bfec39' in str(relpath):
            RGIpath = relpath + '/'
            break

    args = parseCommandLine()
    createDirs()

    roughGenomesList = getFilesList(args.input)
    genomesList = checkFiles(roughGenomesList)
    GENOMES, FILENAMES = clearGlobalDicts()

    resultsDict = getResults(genomesList, RGIpath)
    #print resultsDict

