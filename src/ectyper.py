#!/usr/bin/env python

"""
    Predictive genomics for _E. coli_ from the command line.
    Currently includes serotyping and VF finding.
"""

import logging
import definitions
import src.blastFunctions
import src.commandLineOptions
import src.genomeFunctions
import src.loggingFunctions
import src.resultsToTable
import json
import csv
import os.path
import sys

log_file = src.loggingFunctions.initialize_logging()
log = logging.getLogger(__name__)


def run_program():
    """
    Wrapper for both the serotyping and virulence finder
    The program needs to do the following
    (1) Get names of all genomes being tested
    (2) Create a BLAST database of those genomes
    (3) Query for serotype and/or virulence factors
    (4) Parse the results
    (5) Display the results
    :return: success or failure

    """

    log.info("Starting ectyper. Logfile is: " + str(log_file))
    args = src.commandLineOptions.parse_command_line()
    log.debug(args)

    # run none / one / both of serotyper and vffinder
    query_file = None

    if args.serotyper and args.virulenceFactors:
        query_file = definitions.SEROTYPE_AND_VF_FILE
    elif args.virulenceFactors:
        query_file = definitions.VF_FILE
    elif args.serotyper:
        query_file = definitions.SEROTYPE_FILE
    else:
        log.warning('No analyses selected to run. Exiting. Please select one'
                    ' or both of `--serotyper` and `--virulenceFactors`')
        exit(1)

    log.info("Gathering genome file names")
    raw_genome_files = src.genomeFunctions.get_files_as_list(args.input)
    log.debug(raw_genome_files)

    log.info("Gathering genome names from files")
    (all_genomes_list,
     all_genomes_files) = src.genomeFunctions.get_genome_names_from_files(
        raw_genome_files)
    log.debug(all_genomes_list)
    log.debug(all_genomes_files)

    log.info("Creating blast database")
    blast_db = src.blastFunctions.create_blast_db(all_genomes_files)

    log.info("Blast queries %s against the database of input files",
             query_file)
    #blast_output_file = src.blastFunctions.run_blast(query_file, blast_db)
    blast_output_file = os.path.abspath('C:\\Users\gnial\AppData\Local\Temp\ectyper_blastdb.output')

    log.info("Parsing blast results in %s", blast_output_file)
    # We want to make the parsing function generalizable, not dependent
    # on testing for serotype, vfs etc. specifically
    # The parser will apply any functions given to it to a blast record
    list_of_parsing_dict = src.genomeFunctions.get_list_of_parsing_dict(
        args)

    parsed_results = src.blastFunctions.parse_blast_results(args,
                                                            blast_output_file,
                                                            list_of_parsing_dict)
    #print the requested format
    if args.tabular:
        log.info("Printing results in tabular format")

        writer = csv.writer(sys.stdout, delimiter="\t")
        writer.writerows(src.resultsToTable.results_dict_to_table(all_genomes_files
                                                                  ,all_genomes_list
                                                                  ,parsed_results))
    else:
        log.info("Printing results in JSON format")
        print(json.dumps(parsed_results))

    log.info("Done")
