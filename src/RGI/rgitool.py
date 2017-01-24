#!/usr/bin/env python

import json, sys, os, argparse, logging

TEMP_SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../Serotyper/'))
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../'))

from ecvalidatingfiles import *
from sharedmethods import *

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"

GENOMES = {}
GENOMENAMES = {}

def getGENOMESDict():
    return GENOMES

def setGlobalDicts(genomes_dict, genomenames_dict):
    global GENOMENAMES
    global GENOMES
    GENOMENAMES = genomenames_dict
    GENOMES = genomes_dict


def parseCommandLine():
    """
    Initalizing the main commands of the command line for the project.
    - input: refers to the location of the file(s) that will be processed
    - out: refers to the output of the program. Default is STDOUT
    - tsv: if the user wants a csv copy of the results
    - min: the minimum number of genomes containing a certain gene
    - p: trigger for for filtering strict matches

    :return parser.parse_args(): Data from the commands.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-in", "--input", help="Location of new file(s). Can be a single file or a directory", required=True)
    parser.add_argument("-out", "--output", type=argparse.FileType('w'), help="Output of the program. Default is STDOUT.", default=sys.stdout)
    parser.add_argument("-csv", type=int, help="If set to 1, the results will be sent to a .csv file in the temp/Results folder. Options are 0 and 1, default=1.", default=1, choices=[0,1])
    parser.add_argument("-min", "--minGenomes", type=int, help="Minimum number of genomes threshold for a gene. Default is 1.", default=1)
    parser.add_argument("-p", "--perfectMatches", type=int, help="If set to 1, the result of the AMR tool will contain only perfect matches (no strict). Options are 0 and 1, default is 0.", default=0, choices=[0,1])

    return parser.parse_args()


def getResults(genomesList, RGIpath):
    """
    Going through each file/genome entered and running the McMaster RGI tool to get the results and transfer them in a TSV file.

    :param genomesList: list of files containing genomes
    :param RGIpath: path to the McMaster RGI tool
    """

    logging.info('Entering the CARD RGI tool.')

    global GENOMES
    global GENOMENAMES

    rel_dir = 'temp/Results/RGI/'

    for genome_file in genomesList:
        filename = os.path.basename(genome_file)
        filename = os.path.splitext(filename)
        genome_name = GENOMENAMES[filename[0]]
        genome_name = getGenomeName(genome_name, filename[0])

        logging.info('Getting results for ' + str(filename[0]))

        #Generating the results from the RGI tool
        out = rel_dir + genome_name
        temp_result = subprocess.call(['rgi.py', "-i", genome_file, "-o", out])

        #Obtaining a readable JSON file of the results
        formatted_out = rel_dir + genome_name + '_FORMATTED'
        temp_result = subprocess.call(['rgi_jsonformat', "-i", SCRIPT_DIRECTORY + '../../' + out + '.json', "-o", formatted_out])

        with open(SCRIPT_DIRECTORY + '../../' + formatted_out + '.json', 'r') as temp_file:
            GENOMES[genome_name] = json.load(temp_file)

        #Generating a CSV file of the results
        csv_out = rel_dir + genome_name
        temp_result = subprocess.call(['rgi_jsontab', "-i", SCRIPT_DIRECTORY + '../../' + formatted_out + '.json', "-o", csv_out])

        os.rename(SCRIPT_DIRECTORY + '../../' + csv_out + '.txt', SCRIPT_DIRECTORY + '../../' + csv_out + '.tsv')
        os.remove(SCRIPT_DIRECTORY + '../../' + formatted_out + '.json')
        os.remove(SCRIPT_DIRECTORY + '../../' + out + '.json')

    logging.info('Exiting the McMaster RGI tool.')


def getPerfectMatches(genomesDict):
    """
    Going through the results dictionary to filter out the entries that aren't perfect.

    :param genomesDict: dictionary containing raw results (unfiltered).
    """

    global GENOMES
    for genome_name, genome_info in genomesDict.iteritems():
        GENOMES[genome_name] = {}
        for contig_name, contig_info in genome_info.iteritems():
            GENOMES[genome_name][contig_name] = {}
            for id, id_info in contig_info.iteritems():
                if isinstance(id_info, dict)\
                and 'type_match' in id_info.keys()\
                and id_info['type_match'] == 'Perfect':
                    GENOMES[genome_name][contig_name][id] = id_info



def getGeneDict(genomesDict):
    """
    Collecting only gene names and their presence by going through the results dictionary.

    :param genomesDict: results dictionary
    """
    logging.info('Creating the filtered result dictionary.')
    global GENOMES
    GENOMES = {}

    for genome_name, genome_info in genomesDict.iteritems():
        GENOMES[genome_name] = {}
        for contig_name, contig_info in genome_info.iteritems():
            for id, id_info in contig_info.iteritems():
                if isinstance(id_info, dict) and 'ARO_name' in id_info.keys():
                    gene_name = id_info['ARO_name']
                    GENOMES[genome_name][gene_name] = 1


if __name__ == '__main__':

    logging.basicConfig(filename=SCRIPT_DIRECTORY + 'rgi.log',level=logging.INFO)
    logging.info("This log file doesn't contain the logs from the actual RGI tool. To see these logs, go to your downloaded RGI folder and open the app.log file.")


    #Retrieving the RGI python path.
    RGIpath = ''
    for relpath in sys.path:
        if 'release-rgi-v3.1.1-58cad6a3b443abb290cf3df438fe558bc5bfec39' in str(relpath):
            RGIpath = relpath + '/'
            break


    args = parseCommandLine()

    if args.input == None:
        logging.info('No inputs were given.')
        print 'Error'

    else:
        createDirs()

        roughGenomesList = getFilesList(args.input)
        genomesList = checkFiles(roughGenomesList)

        if isinstance(genomesList, list):
            GENOMES = getGENOMES()
            GENOMENAMES = getGENOMENAMES()

            clearGlobalDicts()

            getResults(genomesList, RGIpath)
            if args.perfectMatches == 1:
                logging.info('Getting perfect matches for ' + str(GENOMES))
                getPerfectMatches(GENOMES)

            getGeneDict(GENOMES)
            resultDict = filterGenes(GENOMES, args.minGenomes)

            if args.csv == 1:
                toTSV(resultDict, 'RGI_Results')

            print resultDict
        else:
            print genomesList
