#!/usr/bin/env python

import argparse
import os
import subprocess
import logging
import re
from ecprediction import getProductPercentage

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO


SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
GENOMES = {}

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

    parser.add_argument("-input", help="Location of new file(s). Can be a single file or a directory")
    parser.add_argument("-out", "-output", help="Output of the program. Default is STDOUT.", default='STDOUT')
    parser.add_argument("-pi", "-percentIdentity", type=int, help="Percentage of identity wanted to use against the database. From 0 to 100, default is 90%.", default=90)
    parser.add_argument("-pl", "-percentLength", type=int, help="Percentage of length wanted to use against the database. From 0 to 100, default is 90%.", default=90)
    parser.add_argument("-v", "-verbose", help="Information desired, true being full information, false being the serotype only.", default="false")
    parser.add_argument("-csv", help="If set to true, the results will be sent to a .csv file in the temp/Results folder.", default='true')

    return parser.parse_args()


def getFilesList(data):
    """
    Creating a list out of the files entered (where each file name is its absolute path). This creates a uniform
    format that works for both single files and directories.

    :param data: Data (file or directory) taken from the input command.
    :return filesList: List of all the files found in the data.
    """
    print data

    filesList = []

    if os.path.isdir(data):
        print "It's a directory"
        logging.info("Using files from " + data)

        for root, dirs, files in os.walk(data):
            for filename in files:
                filesList.append(os.path.join(root,filename))

    else:
        print "It's a file"
        logging.info("Using file " + data)
        filesList.append(os.path.abspath(data))

    return sorted(filesList)


def getGenomeName(recordID, filename):
    """
    Getting the name of the genome by hierarchy to store or search in the GENOMES dictionary in the later called methods.

    :param recordID: ID of a record of a sequence.
    :param filename: Name of the file containing the record.
    :return genomeName: Name of the genome contained in the file (or sequence).
    """

    recordID = str(recordID)

    if re.search('lcl\|([\w-]*)', recordID):
        match = re.search('lcl\|([\w-]*)', recordID)
        match = str(match.group())
        genome_name = match.split('|')[1]

    elif re.search('(^[a-zA-Z][a-zA-Z]\w{6}\.\d)',recordID):
        match = re.search('(\w{8}\.\d)',recordID)
        genome_name = str(match.group())

    elif re.search('(ref\|\w{2}_\w{6}|gb\|\w{8}|emb\|\w{8}|dbj\|\w{8})',recordID):
        match = re.search('(ref\|\w{2}_\w{6}|gb\|\w{8}|emb\|\w{8}|dbj\|\w{8})',recordID)
        match = str(match.group())
        genome_name = match.split('|')[1]

    elif re.search('gi\|\d{8}', recordID):
        match = re.search('gi\|\d{8}', recordID)
        match = str(match.group())
        genome_name = match.split('|')[1]

    else:
        genome_name = filename

    return genome_name


def checkFiles(genomesList):
    """
    Creating a list containing only valid .fasta files out of the previously formed list (from the getListGenomes method).
    If the newly created list is empty, the program exits with a warning.
    This filters out the invalid files.

    :param genomesList: Unfiletered files list.
    :return newGenomesList: List containing only valid fasta files.
    """
    global GENOMES
    newGenomesList = []

    for file in genomesList:
        flag = 0
        for record in SeqIO.parse(file, "fasta"):
            match = re.search('(^[a-zA-Z]+)', str(record.seq))
            if not match:
                flag = 0
                break
            else:
                flag = 1

        if flag>0:
            newGenomesList.append(file)

            filename = os.path.basename(file)
            filename = os.path.splitext(filename)

            for record in SeqIO.parse(file,"fasta"):
                genome_name = getGenomeName(record, filename)
                if not genome_name in GENOMES:
                    GENOMES[genome_name] = ''

        else:
            logging.warning("File " + file + " is in invalid format")

    if not newGenomesList:
        logging.error("No valid fasta files \n Exiting")
        return 'Error'

    else:
        return sorted(newGenomesList)


def initializeDB():
    """
    Generating the database if it was not already generated. The database can be found in the tmp folder.

    :return int 0 or 1: 0 being that the database was created successfully (or already existed).
    """

    REL_DIR = SCRIPT_DIRECTORY + '../temp/database/'

    if not os.path.isdir(REL_DIR):
        os.mkdir(REL_DIR)

    if os.path.isfile(REL_DIR + 'ECTyperDB.nin'):
        return 0
    else:
        return subprocess.call(["/usr/bin/makeblastdb", "-in", SCRIPT_DIRECTORY + "../Data/EcOH.fasta ", "-dbtype", "nucl", "-title", "ECTyperDB", "-out", REL_DIR + "ECTyperDB"])


def runBlastQuery(genomesList, db_name):
    """
    Generating the .xml files containing the genomes that were queried by using the command line.

    :param genomesList: List containing all the files to be compared to the database.
    :param db_name: Name of the database to search.
    :return resultList: List of all the new .xml files containing the results of the search.
    """

    REL_DIR = SCRIPT_DIRECTORY + '../temp/xml/'
    resultsList = []

    if not os.path.exists(REL_DIR):
        os.mkdir(REL_DIR)

        if os.path.isdir(REL_DIR):
            print str(REL_DIR) + " is now a directory."


    for file in genomesList:
        filename = os.path.basename(file)
        filename = os.path.splitext(filename)


        newFilename = os.path.abspath(REL_DIR  + filename[0] + '.xml')

        blastn_cline = NcbiblastnCommandline(cmd="blastn", query=file, db= REL_DIR + '../database/' +  db_name, outfmt=5, out= newFilename)

        stdout, stderr = blastn_cline()
        resultsList.append(newFilename)

    logging.info("Generated " + str(len(resultsList)) + " .xml file(s)")
    return sorted(resultsList)


def parseResults(resultsList):
    """
    Searching the result list to find the hsps necessary to identify the sequences.
    The method stores them in a dictionary that will be stored in GENOMES dictionary.

    :param resultsList: List of all the .xml files to be parsed.
    :return alignmentsDict: Dictionary containing all the information about the matches found in the database for each
    genome.
    """

    global GENOMES
    logging.info("Parsing results from " + str(resultsList))

    for result in resultsList:
        result_handle = open(result)
        blast_records = NCBIXML.parse(result_handle)
        filename = os.path.basename(result)
        filename = os.path.splitext(filename)

        for blast_record in blast_records:

            genome_name = getGenomeName(blast_record.query, filename)
            alignmentsDict = {}
            if genome_name in GENOMES:
             alignmentsDict = dict(GENOMES[genome_name])

            for alignment in blast_record.alignments:
                # hsp_length = abs(1-(1-float(alignment.hsps[0].positives)/alignment.length))
                # hsp_identity = float(alignment.hsps[0].identities)/alignment.hsps[0].align_length
                # print str(alignment.title) + ": " + str(hsp_length) + ", " + str(hsp_identity) + ", " + str(alignment.length) + "\n"
                alignmentsDict[alignment.title] = {alignment.length : alignment.hsps[0]}
                GENOMES[genome_name] = alignmentsDict


    return GENOMES