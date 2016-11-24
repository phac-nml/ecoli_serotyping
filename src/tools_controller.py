#!/usr/bin/env python

import subprocess, sys, os, argparse, logging, csv, pandas as pd

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"

serotyper = SCRIPT_DIRECTORY + '../temp/Results/Serotyper_Results.csv'
virulence_factors = SCRIPT_DIRECTORY + '../temp/Results/VF_Results.tsv'
amr = SCRIPT_DIRECTORY + '../temp/Results/RGI_Results.tsv'


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

    return parser.parse_args()


def createReport(tsv_files):
    """
    Generating a TSV file containing all the results from all the tools.

    :param tsv_files: List of applicable TSV files.
    """

    logging.info('Transfering all the result to file Summary.tsv. This file can be found in temp/Results/ folder.')
    if len(tsv_files) == 1:
        filename = os.path.basename(tsv_files[0])
        os.rename(SCRIPT_DIRECTORY + '../temp/Results/' + filename, SCRIPT_DIRECTORY + '../temp/Results/Summary.tsv')

    elif len(tsv_files) == 2:
         file1 = pd.read_csv(tsv_files[0])
         file2 = pd.read_csv(tsv_files[1])
         merged = file1.merge(file2, on='Genome')
         merged.to_csv(SCRIPT_DIRECTORY + '../temp/Results/Summary.tsv', index=False)
    else:
        file1 = pd.read_csv(tsv_files[0])
        file2 = pd.read_csv(tsv_files[1])
        file3 = pd.read_csv(tsv_files[2])
        merged = file1.merge(file2, on='Genome')
        merged.to_csv(SCRIPT_DIRECTORY + '../temp/Results/temp.tsv', index=False)

        file4 = pd.read_csv(SCRIPT_DIRECTORY + '../temp/Results/temp.tsv')
        merged = file4.merge(file3, on='Genome')
        merged.to_csv(SCRIPT_DIRECTORY + '../temp/Results/Summary.tsv', index=False)
        os.remove(SCRIPT_DIRECTORY + '../temp/Results/temp.tsv')


if __name__=='__main__':

    logging.basicConfig(filename=SCRIPT_DIRECTORY + 'controller.log',level=logging.INFO)

    args = parseCommandLine()
    print args
    logging.info('Starting the controller.')

    tsv_files = []

    if args.serotyper != 1 and args.virulencefactors !=1 and args.amr !=1:
        logging.error('No tools (Serotyper, Virulence Factors tool or AMR (RGI) tool) has been selected. Exiting the program.')
        print 'Error'

    else:
        if args.serotyper == 1:
            logging.info('Triggering the use of the Serotyper tool.')
            tsv_files.append(serotyper)
            temp = subprocess.call([SCRIPT_DIRECTORY + "Serotyper/ectyper.py", "--input", args.input, "-pl", str(args.percentLength),
                                    "-pi", str(args.percentIdentity), '-v', str(args.serotypeverbose), '-csv', '1'])

        if args.virulencefactors == 1:
            logging.info('Triggering the use of the Virulence Factors tool.')
            tsv_files.append(virulence_factors)
            temp = subprocess.call([SCRIPT_DIRECTORY + "Virulence_Factors/virulencefactors.py", "--input", args.input, "-pl", str(args.percentLength),
                                    "-pi", str(args.percentIdentity), '-tsv', '1', '-min', str(args.mingenomes)])

        if args.amr == 1:
            logging.info('Triggering the use of the RGI tool.')
            tsv_files.append(amr)
            temp = subprocess.call([SCRIPT_DIRECTORY + "RGI/rgitool.py", "--input", args.input,
                                    '-tsv', '1', '-min', str(args.mingenomes)])

        createReport(tsv_files)


