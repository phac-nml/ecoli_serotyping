#!/usr/bin/env python

import json
import sys
import os

TEMP_SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../Serotyper/'))
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../'))

from ecvalidatingfiles import *
from formatresults import *
from createdirs import *

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"

GENOMES = {}
GENOMENAMES = {}

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
    parser.add_argument("-tsv", help="If set to 1, the results will be sent to a .tsv file in the temp/Results folder. Options are 0 and 1, default=1.", default=1)

    return parser.parse_args()

def getResults(genomesList, RGIpath):
    global GENOMES
    global GENOMENAMES

    rel_dir = 'temp/Results/RGI/'

    for genome_file in genomesList:
        filename = os.path.basename(genome_file)
        filename = os.path.splitext(filename)
        genome_name = GENOMENAMES[filename[0]]
        genome_name = getGenomeName(genome_name, filename[0])

        out = rel_dir + genome_name
        temp_result = subprocess.call(['python', RGIpath + 'rgi.py', "-i", genome_file, "-o", out])

        formatted_out = rel_dir + genome_name + '_FORMATTED'
        temp_result = subprocess.call(['python',  RGIpath + 'formatJson.py', "-i", SCRIPT_DIRECTORY + '../../' + out + '.json', "-o", formatted_out])

        with open(SCRIPT_DIRECTORY + '../../' + formatted_out + '.json', 'r') as temp_file:
            GENOMES[genome_name] = json.load(temp_file)

        csv_out = rel_dir + genome_name
        temp_result = subprocess.call(['python',  RGIpath + 'convertJsonToTSV.py', "-i", SCRIPT_DIRECTORY + '../../' + formatted_out + '.json', "-o", csv_out])
        os.rename(SCRIPT_DIRECTORY + '../../' + csv_out + '.txt', SCRIPT_DIRECTORY + '../../' + csv_out + '.tsv')

        os.remove(SCRIPT_DIRECTORY + '../../' + formatted_out + '.json')
        os.remove(SCRIPT_DIRECTORY + '../../' + out + '.json')



def filterResults(genomesDict):
    global GENOMES
    GENOMES = {}

    for genome_name, genome_info in genomesDict.iteritems():
        GENOMES[genome_name] = {}
        for contig_name, contig_info in genome_info.iteritems():
            for id, id_info in contig_info.iteritems():
                if isinstance(id_info, dict) and 'ARO_name' in id_info.keys():
                    gene_name = id_info['ARO_name']
                    GENOMES[genome_name][gene_name] = 1

def getGENOMES():
    return GENOMES

def setGENOMES(genomes_dict, genome_filenames):
    global GENOMES
    global GENOMENAMES
    GENOMES = genomes_dict
    GENOMENAMES = genome_filenames


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
    GENOMES, useless_dict, GENOMENAMES = clearGlobalDicts()

    getResults(genomesList, RGIpath)
    filterResults(GENOMES)

    if args.tsv == 1:
        toTSV(GENOMES, 'RGI_Results')
    print GENOMES

