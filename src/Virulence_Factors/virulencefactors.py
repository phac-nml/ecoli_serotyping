#!/usr/bin/env python
from __future__ import print_function, division
import os
import sys
import logging
import argparse
import subprocess
import shutil
import json

TEMP_SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../Serotyper/'))
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../'))

from sharedmethods import *
from ecvalidatingfiles import *

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
GENOMES = {}
FILENAMES = {}


def parseCommandLine():
    """
    Initalizing the main commands of the command line for the project.
    - input: refers to the location of the file(s) that will be processed
    - out: refers to the output of the program. Default is STDOUT
    - pi: refers to the percentage of identity wanted. Default is 90%
    - pl: refers to the percentage of length wanted. Default is 90%.
    - csv: if the user wants a csv copy of the results
    - min: the minimum number of genomes containing a certain gene

    :return parser.parse_args(): Data from the commands.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-in", "--input", help="Location of new file(s). Can be a single file or a directory", required=True)
    parser.add_argument("-out", "--output", type=argparse.FileType('w'),
                        help="Output of the program. Default is STDOUT.", default=sys.stdout)
    parser.add_argument("-pi", "--percentIdentity", type=int,
                        help="Percentage of identity wanted to use against the database. From 0 to 100, default is 90%.", default=90, choices=range(0, 100))
    parser.add_argument("-pl", "--percentLength", type=int,
                        help="Percentage of length wanted to use against the database. From 0 to 100, default is 90%.", default=90, choices=range(0, 100))
    parser.add_argument(
        "-csv", type=int, help="If set to 1, the results will be sent to a .csv file in the temp/Results folder. Options are 0 and 1, default is 1.", default=1, choices=[0, 1])
    parser.add_argument("-min", "--minGenomes", type=int,
                        help="Minimum number of genomes threshold for a virulence factor. Default is 1.", default=0)

    return parser.parse_args()


def initializeDB():
    """
    Generating the database if it was not already generated. The database can be found in the tmp folder.

    :return int 0 or 1: 0 being that the database was created successfully (or already existed).
    """

    REL_DIR = os.path.join(
        SCRIPT_DIRECTORY, '../../temp/databases/VF_Database/')

    db_in = os.path.join(SCRIPT_DIRECTORY, '../../Data/repaired_ecoli_vfs.ffn')
    db_out = os.path.join(REL_DIR, 'VirulenceFactorsDB')

    makeblastdb = ("makeblastdb", "-in", db_in, "-dbtype", "nucl",
                   "-title", "VirulenceFactorsDB", "-out", db_out)

    if os.path.isfile(REL_DIR + 'VirulenceFactorsDB.nin'):
        logging.info('Database already exists.')
        retcode = 0

    else:
        logging.info('Generating the database.')

        retcode = subprocess.call(makeblastdb)

    return retcode


def searchDB(genomesList):
    """
    Querying the database with a single file combining all the genomes.

    :param genomesList: files that will be used against the database
    :return new_filename: .xml file containing the results from querying the database.
    """

    REL_DIR = os.path.join(SCRIPT_DIRECTORY, '../../temp/xml/')

    if len(genomesList) > 1:
        combined_genomes = os.path.join(SCRIPT_DIRECTORY,
                                        '../../temp/Uploads/combined_genomesVF.fasta')

        # Copying the content of every fasta file into one file to simplify the
        # database search
        with open(combined_genomes, 'wb') as outfile:
            for file in genomesList:
                with open(file, 'rb') as fastafile:
                    shutil.copyfileobj(fastafile, outfile, 1024 * 1024 * 10)

        new_filename = os.path.abspath(
            os.path.join(REL_DIR, 'combined_genomesVF.xml'))

    else:
        filename = os.path.basename(genomesList[0])
        filename, ext = os.path.splitext(filename)
        combined_genomes = genomesList[0]
        new_filename = os.path.abspath(
            os.path.join(REL_DIR, filename + '.xml'))

    # Querying the database

    db_path = os.path.join(
        REL_DIR, '../databases/VF_Database/VirulenceFactorsDB')

    blastn_cline = NcbiblastnCommandline(cmd="blastn", query=combined_genomes,
                                         db=db_path, outfmt=5, out=new_filename)
    stdout, stderr = blastn_cline()

    logging.info("Searched the database.")
    return new_filename


def parseFile(result_file, perc_len, perc_id):
    """
    Going through each record in the result file and filter them by the percentages provided.

    :param result_file: file containing all the results from querying the database.
    :param perc_len: percentage length cutoff
    :param perc_id: percentage identity cutoff
    """

    global GENOMES
    global FILENAMES
    logging.info("Parsing results from " + str(result_file))

    result_handle = open(result_file)
    blast_records = NCBIXML.parse(result_handle)

    perc_len = perc_len / 100
    perc_id = perc_id / 100

    for blast_record in blast_records:
        filename = FILENAMES[blast_record.query]
        genome_name = getGenomeName(blast_record.query, filename)
        alignmentsDict = {}
        if genome_name in GENOMES:
            alignmentsDict = dict(GENOMES[genome_name])

        # Find the gene name + compare with the identity and length cutoffs.
        for alignment in blast_record.alignments:
            match = re.search('(\s\(.*\)\s)', alignment.title)
            match = match.group()
            match = match.split('(')[1]
            match = match.split(')')[0]
            align_title = str(match)

            hsp = alignment.hsps[0]

            tmp_perc_len = abs(1 - (1 - hsp.positives) / alignment.length)

            tmp_perc_id = hsp.positives / hsp.align_length

            # standard check for perc_len and perc_id
            if tmp_perc_len > perc_len and tmp_perc_id > perc_id:

                # we do this after the perc_len/id checks so as to not incl
                # blank dicts if they arent valid
                contig_accession = blast_record.query.split(' ')[0]
                # check if this contig already has a gene assoc with it
                if contig_accession not in alignmentsDict.keys():
                    alignmentsDict[contig_accession] = []

                # intend to make it: {filename{contig_accession[{'START','STOP','ORIENTATION','GENE_NAME'}]}}
                # we use a list after contig accession just in case that contig
                # has more than 1 genes assoc with it
                geneDict = {}
                geneDict['GENE_NAME'] = align_title
                # check if the subject (what we're looking for) was found in forward or reverse
                # this means found in forward
                if hsp.sbjct_start < hsp.sbjct_end:
                    geneDict['ORIENTATION'] = '+'
                else:
                    geneDict['ORIENTATION'] = '-'

                geneDict['START'] = hsp.query_start
                geneDict['STOP'] = hsp.query_end

                # dupl check before adding it in, we do a full comparison in case the same gene is found in two different locations within a contig
                if not any(h == geneDict for h in alignmentsDict[contig_accession]):
                    alignmentsDict[contig_accession].append(geneDict)

            GENOMES[genome_name] = alignmentsDict


if __name__ == '__main__':

    log_path = os.path.join(SCRIPT_DIRECTORY, 'virulencefactors.log')
    logging.basicConfig(filename=log_path, level=logging.INFO)

    args = parseCommandLine()
    createDirs()

    logging.info('Starting Virulence Factors tool')

    roughGenomesList = getFilesList(args.input)
    genomesList = checkFiles(roughGenomesList)

    GENOMES = getGENOMES()
    FILENAMES = getFILENAMES()
    clearGlobalDicts()

    if isinstance(genomesList, list):

        if initializeDB() == 0:

            results_file = searchDB(genomesList)
            parseFile(results_file, args.percentLength, args.percentIdentity)
            #resultsDict = filterGenes(GENOMES, args.minGenomes)
            resultsDict = GENOMES

            if args.csv == 1:
                toTSV(resultsDict, 'VF_Results')

            json.dump(resultsDict, sys.stdout)
            logging.info('Program ended successfully.')

    else:
        print(genomesList)
