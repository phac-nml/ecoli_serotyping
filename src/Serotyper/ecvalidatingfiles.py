#!/usr/bin/env python

import argparse
import os
import subprocess
import logging
import sys
import shutil

TEMP_SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../'))

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from sharedmethods import *

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
GENOMES = {}
FILENAMES = {}



def setGlobalDicts():
    global GENOMES
    global GENOMENAMES
    global FILENAMES
    GENOMES = getGENOMES()
    FILENAMES = getFILENAMES()

def setFILENAMESDict(filenames_dict):
    global FILENAMES
    FILENAMES = filenames_dict


def parseCommandLine():
    """
    Initalizing the main commands of the command line for the project.
    - input: refers to the location of the file(s) that will be processed
    - out: refers to the output of the program. Default is STDOUT
    - pi: refers to the percentage of identity wanted. Default is 90%
    - pl: refers to the percentage of length wanted. Default is 90%.
    - v: refers to the verbosity.
    - csv: if the user wants a csv copy of the results.

    :return parser.parse_args(): Data from the commands.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-in", "--input", help="Location of new file(s). Can be a single file or a directory", required=True)
    parser.add_argument("-out", "--output", type=argparse.FileType('w'), help="Output of the program. Default is STDOUT.", default=sys.stdout)
    parser.add_argument("-pi", "--percentIdentity", type=int, help="Percentage of identity wanted to use against the database. From 0 to 100, default is 90%.", default=90, choices=range(0,100))
    parser.add_argument("-pl", "--percentLength", type=int, help="Percentage of length wanted to use against the database. From 0 to 100, default is 90%.", default=90, choices=range(0,100))
    parser.add_argument("-v", "--verbose", type=int, help="Information desired, 1 being full information, 0 being the serotype only. Options are 0 and 1, default is 0", default=0, choices=[0,1])
    parser.add_argument("-csv", type=int, help="If set to 1, the results will be sent to a .csv file in the temp/Results folder. Options are 0 and 1, default is 1.", default=1, choices=[0,1])
    parser.add_argument("-g", "--galaxy", type=int, help="Galaxy is in use", default=0, choices=[0,1])
    parser.add_argument("-outfmt", "--outputFormat", help="Format of the results (Table or JSON).", default="json")

    return parser.parse_args()


def initializeDB():
    """
    Generating the database if it was not already generated. The database can be found in the tmp folder.

    :return int 0 or 1: 0 being that the database was created successfully (or already existed).
    """

    REL_DIR = SCRIPT_DIRECTORY + '../../temp/databases/Serotyping_Database/'

    if os.path.isfile(REL_DIR + 'ECTyperDB.nin'):
        logging.info('Database exists.')
        return 0
    else:
        logging.info('Generating the database.')
        return subprocess.call(["/usr/bin/makeblastdb", "-in", SCRIPT_DIRECTORY + "../../Data/EcOH.fasta ", "-dbtype", "nucl", "-title", "ECTyperDB", "-out", REL_DIR + "ECTyperDB"])


def runBlastQuery(genomesList, db_name):
    """
    Generating the .xml files containing the genomes that were queried by using the command line.

    :param genomesList: List containing all the files to be compared to the database.
    :param db_name: Name of the database to search.
    :return new_filename: .xml file containing the results of the search.
    """

    REL_DIR = SCRIPT_DIRECTORY + '../../temp/xml/'

    if len(genomesList) >1:
        combined_genomes = SCRIPT_DIRECTORY + '../../temp/Uploads/combined_genomes.fasta'

        #Copy each file content into one file to simplify the database search
        with open(combined_genomes, 'wb') as outfile:
            for file in genomesList:
                with open(file, 'rb') as fastafile:
                    shutil.copyfileobj(fastafile, outfile,1024*1024*10)

        new_filename = os.path.abspath(REL_DIR  + 'combined_genomes.xml')

    else:
        filename = os.path.basename(genomesList[0])
        filename = os.path.splitext(filename)
        combined_genomes = genomesList[0]
        new_filename = os.path.abspath(REL_DIR + str(filename[0]) + '.xml')

    #Run the BLAST query
    blastn_cline = NcbiblastnCommandline(cmd="blastn", query=combined_genomes, db= REL_DIR + '../databases/Serotyping_Database/' +  db_name, outfmt=5, out= new_filename)
    stdout, stderr = blastn_cline()

    logging.info("Searched the database.")
    return new_filename


def parseResults(result_file):
    """
    Searching the result list to find the hsps necessary to identify the sequences.
    The method stores them in a dictionary that will be stored in GENOMES dictionary.

    :param resultsList: List of all the .xml files to be parsed.
    :return alignmentsDict: Dictionary containing all the information about the matches found in the database for each
    genome.
    """

    global GENOMES
    global FILENAMES
    logging.info("Parsing results from " + str(result_file))

    result_handle = open(result_file)
    blast_records = NCBIXML.parse(result_handle)

    for blast_record in blast_records:
        #Obtain record information
        filename = FILENAMES[blast_record.query]
        genome_name = getGenomeName(blast_record.query, filename)

        alignmentsDict = {}
        if genome_name in GENOMES:
            alignmentsDict = dict(GENOMES[genome_name])

        for alignment in blast_record.alignments:
            alignmentsDict[alignment.title] = {alignment.length : alignment.hsps[0]}
            GENOMES[genome_name] = alignmentsDict

    return GENOMES