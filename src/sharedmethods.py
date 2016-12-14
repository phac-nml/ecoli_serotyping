#!/usr/bin/env python

import os
import logging
import re
import csv

from Bio import SeqIO
SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"

GENOMES = {}
FILENAMES = {}
GENOMENAMES = {}

def getGENOMES():
    return GENOMES

def getFILENAMES():
    return FILENAMES

def getGENOMENAMES():
    return GENOMENAMES

def setGENOMES(genomes_dict):
    global GENOMES
    GENOMES = genomes_dict

def setGENOMENAMES(genomenames_dict):
    global GENOMENAMES
    GENOMENAMES = genomenames_dict

def setFILENAMES(filenames_dict):
    global FILENAMES
    FILENAMES = filenames_dict

def clearGlobalDicts():
    global GENOMES
    global FILENAMES
    global GENOMENAMES
    GENOMES = {}
    FILENAMES = {}
    GENOMENAMES= {}


def createDirs():

    if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/'):
        os.mkdir(SCRIPT_DIRECTORY + '../temp/')

    if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/databases/'):
        os.mkdir(SCRIPT_DIRECTORY + '../temp/databases/')

    if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/databases/Serotyping_Database/'):
        os.mkdir(SCRIPT_DIRECTORY + '../temp/databases/Serotyping_Database/')

    if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/databases/VF_Database/'):
        os.mkdir(SCRIPT_DIRECTORY + '../temp/databases/VF_Database/')

    if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/Results/'):
        os.mkdir(SCRIPT_DIRECTORY + '../temp/Results/')

    if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/xml/'):
        os.mkdir(SCRIPT_DIRECTORY + '../temp/xml/')

    if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/Uploads/'):
        os.mkdir(SCRIPT_DIRECTORY + '../temp/Uploads/')

    if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/Results/RGI/'):
        os.mkdir(SCRIPT_DIRECTORY + '../temp/Results/RGI/')


def getFilesList(data):
    """
    Creating a list out of the files entered (where each file name is its absolute path). This creates a uniform
    format that works for both single files and directories.

    :param data: Data (file or directory) taken from the input command.
    :return filesList: List of all the files found in the data.
    """

    filesList = []

    if os.path.isdir(data):
        logging.info("Using files from " + data)

        for root, dirs, files in os.walk(data):
            for filename in files:
                filesList.append(os.path.join(root,filename))

    else:
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
    global FILENAMES
    global GENOMENAMES
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

    if recordID not in FILENAMES:
        FILENAMES[recordID] = filename
    if filename not in GENOMENAMES:
        GENOMENAMES[filename] = recordID

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

    for genome_file in genomesList:
        flag = 0
        for record in SeqIO.parse(genome_file, "fasta"):
            match = re.search('(^[a-zA-Z]+)', str(record.seq))
            if not match:
                flag = 0
                break
            else:
                flag = 1

        if flag>0:
            newGenomesList.append(genome_file)

            filename = os.path.basename(genome_file)
            filename = os.path.splitext(filename)

            for record in SeqIO.parse(genome_file,"fasta"):
                genome_name = getGenomeName(record.description, filename[0])
                if not genome_name in GENOMES:
                    GENOMES[genome_name] = ''

        else:
            logging.warning("File " + genome_file + " is in invalid format")

    if not newGenomesList:
        logging.error("No valid fasta files \n Exiting")
        return 'Error'

    else:
        return sorted(newGenomesList)


def filterGenes(genomesDict, threshold):
    """
    Filtering through the final result dictionary to only keep the genes that meet the threshold requirement.

    :param genomesDict: dictionary containing the final results for each genome
    :param threshold: minimum of presence necessary for a gene to keep it in the dictionary
    :return resultDict: the filtered dictionary
    """

    resultDict = {}
    threshDict = {}

    if threshold == 0:
        return genomesDict

    logging.info('Filtering the result dictionaries with following threshold: ' + str(threshold) + '.')
    for genome_name, gene_info in genomesDict.iteritems():
        if isinstance(gene_info, dict):
            for gene_name in gene_info.keys():
                if gene_name in threshDict and gene_info[gene_name] == 1:
                    threshDict[gene_name] += 1
                elif gene_info[gene_name] == 1:
                    threshDict[gene_name] = 1

    for genome_name, gene_info in genomesDict.iteritems():
        resultDict[genome_name] = {}
        if isinstance(gene_info, dict):
            tempDict = {}
            for gene_name in threshDict.keys():
                if gene_name in gene_info and threshDict[gene_name]>=threshold:
                    tempDict[gene_name] = gene_info[gene_name]
            resultDict[genome_name] = tempDict

    return resultDict


def toTSV(data, results_filename):
    """
    Writing a TSV file containing the result dictionary.

    :param data: result dictionary
    :param results_filename: name of the result file (Virulence_Factors or RGI)
    """
    logging.info('Writing TSV file ' + results_filename + '.tsv. You can find it in the temp/Results/ folder.')

    headers = ['Genome']

    for genome_name, gene_info in data.iteritems():
        if isinstance(gene_info, dict):
            for gene_name in gene_info.keys():
                if gene_name not in headers:
                    headers.append(gene_name)
        else:
            data[genome_name] = {}

    with open(SCRIPT_DIRECTORY + '../temp/Results/' + results_filename +  '.csv', 'wb') as tsvfile:
        tsvwriter = csv.DictWriter(tsvfile, headers)
        tsvwriter.writeheader()

        for genome_name, gene_info in data.iteritems():
            row = {'Genome': genome_name}
            for gene_name in headers[1:]:
                if gene_name in gene_info:
                    row.update({gene_name : gene_info[gene_name]})
                else:
                    row.update({gene_name : '0'})
            tsvwriter.writerow(row)