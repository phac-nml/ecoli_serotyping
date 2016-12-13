#!/usr/bin/env python

import subprocess, sys, os, argparse, logging, csv, pandas as pd, ast

TEMP_SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../Serotyper/'))

from ectyper_formatting import *


SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"

serotyper = SCRIPT_DIRECTORY + '../../temp/Results/Serotyper_Results.csv'
virulence_factors = SCRIPT_DIRECTORY + '../../temp/Results/VF_Results.csv'
amr = SCRIPT_DIRECTORY + '../../temp/Results/RGI_Results.csv'


def parseCommandLine():
    """
    Initalizing the main commands of the command line for the project.
    - input: refers to the location of the file(s) that will be processed
    - s: trigger for the Serotyper
    - vf: trigger for the Virulence Factors tool
    - amr: trigger for the AMR(RGI) tool
    - pi: refers to the percentage of identity wanted. Default is 90%
    - pl: refers to the percentage of length wanted. Default is 90%
    - sv: refers to the verbosity for the Serotyper
    - min: the minimum number of genomes containing a certain gene
    :return parser.parse_args(): Data from the commands.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-in", "--input", help="Location of new file(s). Can be a single file or a directory", required=True)
    parser.add_argument("-s", "--serotyper", type=int, help="Trigger the use of the E. coli Serotyper. Options are 0 and 1. Default is 0 (false).", default=0, choices=[0,1])
    parser.add_argument("-vf", "--virulencefactors",type=int, help="Trigger the use of the Virulence Factors tool. Options are 0 and 1. Default is 0 (false).", default=0, choices=[0,1])
    parser.add_argument("-amr", type=int, help="Trigger the use of the AMR (RGI) tool. Options are 0 and 1. Default is 0 (false).", default=0, choices=[0,1])
    parser.add_argument("-pi", "--percentIdentity", type=int, help="Percentage of identity wanted to use against the database. From 0 to 100, default is 90%.", default=90, choices=range(0,100))
    parser.add_argument("-pl", "--percentLength", type=int, help="Percentage of length wanted to use against the database. From 0 to 100, default is 90%.", default=90, choices=range(0,100))
    parser.add_argument("-sv", "--serotypeverbose", type=int, help="Information desired, 1 being full information, 0 being the serotype only. Options are 0 and 1, default is 0", default=0, choices=[0,1])
    parser.add_argument("-min", "--mingenomes",type=int, help="Minimum number of genomes threshold for a virulence factor. Default is 1.", default=1)
    parser.add_argument("-csv", type=int, help="If set to 1, the results will be sent to a .csv file in the temp/Results folder. Options are 0 and 1, default=1.", default=1, choices=[0,1])
    parser.add_argument("-p", "--perfectMatches", type=int, help="If set to 1, the result of the AMR tool will contain only perfect matches (no strict). Options are 0 and 1, default is 0.", default=0, choices=[0,1])
    parser.add_argument("-avf", "--allVFs", type=int, help="If set to 1, the result of the Virulence Factors tool will contain all virulence factors, not only the top ten. Options are 0 and 1, default is 1.", default=1, choices=[0,1])

    return parser.parse_args()


def createReport(csv_files):
    """
    Generating a TSV file containing all the results from all the tools.

    :param csv_files: List of applicable TSV files.
    """

    logging.info('Transfering all the result to file Summary.csv. This file can be found in temp/Results/ folder.')
    if len(csv_files) == 1:
        filename = os.path.basename(csv_files[0])
        os.rename(SCRIPT_DIRECTORY + '../../temp/Results/' + filename, SCRIPT_DIRECTORY + '../../temp/Results/Summary.csv')

    elif len(csv_files) == 2:
         file1 = pd.read_csv(csv_files[0])
         file2 = pd.read_csv(csv_files[1])
         merged = file1.merge(file2, on='Genome')
         merged.to_csv(SCRIPT_DIRECTORY + '../../temp/Results/Summary.csv', index=False)
    else:
        file1 = pd.read_csv(csv_files[0])
        file2 = pd.read_csv(csv_files[1])
        file3 = pd.read_csv(csv_files[2])
        merged = file1.merge(file2, on='Genome')
        merged.to_csv(SCRIPT_DIRECTORY + '../../temp/Results/temp.csv', index=False)

        file4 = pd.read_csv(SCRIPT_DIRECTORY + '../../temp/Results/temp.csv')
        merged = file4.merge(file3, on='Genome')
        merged.to_csv(SCRIPT_DIRECTORY + '../../temp/Results/Summary.csv', index=False)
        os.remove(SCRIPT_DIRECTORY + '../../temp/Results/temp.csv')


def mergeResults(serotyper_dict, vf_dict, amr_dict):
    """
    Going through all the dictionaries to merge them into a single one.


    :param serotyper_dict: Dictionary containing the Serotyper results.
    :param vf_dict: Dictionary containing the Virulence Factors results.
    :param amr_dict: Dictionary containing the AMR results.
    :return resultDict: Dictionary containing merged results.
    """
    resultsDict = {}

    if bool(vf_dict):
        for genome_name, vf_info in vf_dict.iteritems():
            if genome_name not in resultsDict.keys():
                resultsDict[genome_name] = {}
            resultsDict[genome_name]['Virulence Factors'] = vf_info

    if bool(amr_dict):
        for genome_name, amr_info in amr_dict.iteritems():
            if genome_name not in resultsDict.keys():
                resultsDict[genome_name] = {}
            resultsDict[genome_name]['Antimicrobial Resistance'] = amr_info

    if bool(serotyper_dict):
        for genome_name, serotype_info in serotyper_dict.iteritems():
            if genome_name not in resultsDict.keys():
                resultsDict[genome_name] = {}

            resultsDict[genome_name]['Serotype']= serotype_info

    return resultsDict


if __name__=='__main__':

    logging.basicConfig(filename=SCRIPT_DIRECTORY + 'controller.log',level=logging.INFO)

    args = parseCommandLine()
    logging.info('Starting the controller.')

    csv_files = []
    serotyper_out = '{}'
    vf_out = '{}'
    amr_out = '{}'

    if args.serotyper != 1 and args.virulencefactors !=1 and args.amr !=1:
        logging.error('No tools (Serotyper, Virulence Factors tool or AMR (RGI) tool) has been selected. Exiting the program.')
        print 'Error'

    else:
        if args.serotyper == 1:
            logging.info('Triggering the use of the Serotyper tool.')
            csv_files.append(serotyper)
            serotyper_out = subprocess.check_output([SCRIPT_DIRECTORY + "../Serotyper/ectyper.py", "--input", args.input, "-pl", str(args.percentLength),
                                    "-pi", str(args.percentIdentity), '-v', str(args.serotypeverbose), '-csv', str(args.csv)])

        if args.virulencefactors == 1:
            logging.info('Triggering the use of the Virulence Factors tool.')
            csv_files.append(virulence_factors)
            vf_out = subprocess.check_output([SCRIPT_DIRECTORY + "../Virulence_Factors/virulencefactors.py", "--input", args.input, "-pl", str(args.percentLength),
                                    "-pi", str(args.percentIdentity), '-csv', str(args.csv) , '-min', str(args.mingenomes)])

        if args.amr == 1:
            logging.info('Triggering the use of the RGI tool.')
            csv_files.append(amr)
            amr_out = subprocess.check_output([SCRIPT_DIRECTORY + "../RGI/rgitool.py", "--input", args.input,
                                    '-csv', str(args.csv), '-min', str(args.mingenomes), "-p", str(args.perfectMatches)])
            amr_out = amr_out.split('\n')[-2]

        if args.csv == 1:
            createReport(csv_files)

        resultDict = mergeResults(ast.literal_eval(serotyper_out), ast.literal_eval(vf_out), ast.literal_eval(amr_out))

        print resultDict


