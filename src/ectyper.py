#!/usr/bin/env python

"""
    Predictive genomics for _E. coli_ from the command line.
    Currently includes serotyping and VF finding.
"""

import logging.config
import definitions
import src.blastFunctions
import src.commandLineOptions
import src.genomeFunctions

# need to set disable existing loggers to False, otherwise the module
# logging will not work as intended
logging.config.fileConfig(definitions.LOGGER_CONFIG,
                          disable_existing_loggers=False)
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

    log.info("Starting ectyper")
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
    genome_files = src.genomeFunctions.get_files_as_list(args.input)
    log.debug(genome_files)

    log.info("Gathering genome names from files")
    (all_genomes_list, all_genomes_files) = src.genomeFunctions.get_genome_names_from_files(genome_files)
    log.debug(all_genomes_list)
    log.debug(all_genomes_files)


    log.info("Creating blast database")
    blast_db = src.blastFunctions.create_blast_db(genome_files)

    log.info("Blast queries %s against the database of input files",
             query_file)
    #blast_output_file = src.genomeFunctions.run_blast(query_file, blast_db)
    blast_output_file = '/tmp/tmp7ojkee7v/ectyper_blastdb.output'

    log.info("Parsing blast results in %s", blast_output_file)
    # We want to make the parsing function generalizable, not dependent
    # on testing for serotype, vfs etc. specifically
    # The parser will apply any functions given to it to a blast record
    list_of_parsing_functions = src.genomeFunctions.get_list_of_parsing_functions(args)

    parsed_results = src.blastFunctions.parse_blast_results(args, blast_output_file, list_of_parsing_functions)
    log.info(parsed_results)
    log.info("Done")
