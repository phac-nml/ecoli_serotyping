#!/usr/bin/env python

import os
import sys

TEMP_SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../Serotyper/'))
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../'))

from createdirs import createDirs
from ecvalidatingfiles import *
from formatresults import toTSV

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
GENOMES = {}
FILENAMES = {}

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
    parser.add_argument("-tsv", help="If set to 1, the results will be sent to a .tsv file in the temp/Results folder. Options are 0 and 1, default=1.", default=1)
    parser.add_argument("-min", "--minGenomes", type=int, help="Minimum number of genomes threshold for a virulence factor", default=0)

    return parser.parse_args()


def initializeDB():
    """
    Generating the database if it was not already generated. The database can be found in the tmp folder.

    :return int 0 or 1: 0 being that the database was created successfully (or already existed).
    """

    REL_DIR = SCRIPT_DIRECTORY + '../../temp/databases/VF_Database/'

    if os.path.isfile(REL_DIR + 'VirulenceFactorsDB.nin'):
        return 0
    else:
        return subprocess.call(["/usr/bin/makeblastdb", "-in", SCRIPT_DIRECTORY + "../../Data/repaired_ecoli_vfs.ffn ", "-dbtype", "nucl", "-title", "VirulenceFactorsDB", "-out", REL_DIR + "VirulenceFactorsDB"])


def searchDB(genomesList):

    REL_DIR = SCRIPT_DIRECTORY + '../../temp/xml/'

    if len(genomesList) >1:
        combined_genomes = SCRIPT_DIRECTORY + '../../temp/Uploads/combined_genomesVF.fasta'

        with open(combined_genomes, 'wb') as outfile:
            for file in genomesList:
                with open(file, 'rb') as fastafile:
                    shutil.copyfileobj(fastafile, outfile,1024*1024*10)

        new_filename = os.path.abspath(REL_DIR  + 'combined_genomesVF.xml')

    else:
        filename = os.path.basename(genomesList[0])
        filename = os.path.splitext(filename)
        combined_genomes = genomesList[0]
        new_filename = os.path.abspath(REL_DIR + str(filename[0]) + '.xml')

    blastn_cline = NcbiblastnCommandline(cmd="blastn", query=combined_genomes, db= REL_DIR + '../databases/VF_Database/VirulenceFactorsDB', outfmt=5, out= new_filename)
    stdout, stderr = blastn_cline()

    logging.info("Searched the database.")
    return new_filename


def parseFile(result_file, perc_len, perc_id):
    global GENOMES
    global FILENAMES
    logging.info("Parsing results from " + str(result_file))

    result_handle = open(result_file)
    blast_records = NCBIXML.parse(result_handle)
    perc_len = float(perc_len)/100
    perc_id = float(perc_id)/100

    for blast_record in blast_records:
        filename = FILENAMES[blast_record.query]
        genome_name = getGenomeName(blast_record.query, filename)
        alignmentsDict = {}
        if genome_name in GENOMES:
            alignmentsDict = dict(GENOMES[genome_name])

        for alignment in blast_record.alignments:
            match = re.search('(\s\(.*\)\s)', alignment.title)
            match = match.group()
            match = match.split('(')[1]
            match = match.split(')')[0]
            align_title = str(match)

            tmp_perc_len = abs(1-(1-float(alignment.hsps[0].positives))/alignment.length)
            tmp_perc_id = float(alignment.hsps[0].positives)/alignment.hsps[0].align_length

            if tmp_perc_len > perc_len and tmp_perc_id > perc_id and align_title not in alignmentsDict.keys():
                alignmentsDict[align_title] = 1
                logging.info("")
            else:
                if align_title not in alignmentsDict.keys() or alignmentsDict[align_title]!=1:
                    alignmentsDict[align_title] = 0

            GENOMES[genome_name] = alignmentsDict

def filterVFs(genomesDict, threshold):
    resultDict = {}

    threshDict = {}
    print threshold

    if threshold == 0 or threshold == 1:
        return genomesDict

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
            print str(genome_name)  + '\n---------'
            for gene_name in threshDict.keys():
                if gene_name in gene_info and threshDict[gene_name]>=threshold:
                    tempDict[gene_name] = gene_info[gene_name]
                    print gene_name
            print '_______________'
            resultDict[genome_name] = tempDict

    return resultDict



if __name__=='__main__':

    logging.basicConfig(filename=SCRIPT_DIRECTORY + 'virulencefactors.log',level=logging.INFO)

    args = parseCommandLine()
    createDirs()

    roughGenomesList = getFilesList(args.input)
    genomesList = checkFiles(roughGenomesList)
    GENOMES, FILENAMES, useless_dict = clearGlobalDicts()

    if isinstance(genomesList, list):
        if initializeDB() == 0:
            results_file = searchDB(genomesList)
            parseFile(results_file, args.percentLength, args.percentIdentity)
            resultsDict = filterVFs(GENOMES, args.minGenomes)
            if args.tsv == 1:
                toTSV(resultsDict, 'VF_Results')

            #print resultsDict