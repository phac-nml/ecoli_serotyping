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

    serotype_parsed_results = None
    if args.serotyper:
        serotype_output_file = \
            src.blastFunctions.run_blast(definitions.SEROTYPE_FILE, blast_db)
        serotype_parsed_results = \
            src.blastFunctions.parse_blast_results(
                args,
                serotype_output_file,
                src.genomeFunctions.get_parsing_dict('serotype'))

    virulence_parsed_results = None
    if args.virulenceFactors:
        virulence_output_file = \
            src.blastFunctions.run_blast(definitions.VF_FILE, blast_db)
        virulence_parsed_results = \
            src.blastFunctions.parse_blast_results(
                args,
                virulence_output_file,
                src.genomeFunctions.get_parsing_dict('vf'))

    #parsed_results = None
    if serotype_parsed_results and virulence_parsed_results:
        parsed_results = serotype_parsed_results

        for genome in virulence_parsed_results:
            log.debug("genome: {}".format(genome))
            log.debug("original: {}".format(virulence_parsed_results[
                                                genome]['vf']))
            parsed_results[genome]['vf'] = virulence_parsed_results[genome][
                'vf']
    elif serotype_parsed_results:
        parsed_results = serotype_parsed_results
    else:
        parsed_results = virulence_parsed_results

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
