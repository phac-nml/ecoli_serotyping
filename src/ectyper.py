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
import src.shortReads
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

    # use serotype sequence from `More Serotype`
    query_file = definitions.SEROTYPE_FILE

    raw_genome_files = None

    if args.reads:
        log.info("Gathering genome file names")
        raw_reads_files = src.genomeFunctions.get_files_as_list(args.input)
        log.debug(raw_reads_files)
        _raw_genome_files = []
        # the line below only need to be run once to generate index file for serotyped allele
        # src.shortReads.build_bowtie_index(definitions.SEROTYPE_FILE)
        for raw_reads_file in raw_reads_files:
            # create a fasta file for each of them
            sam_file = src.shortReads.align_reads(raw_reads_file)
            fasta_file = src.shortReads.samToFasta(sam_file)
            _raw_genome_files.append(fasta_file)
        raw_genome_files = _raw_genome_files
    else:
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
    
    serotype_output_file = \
        src.blastFunctions.run_blast(query_file, blast_db)
    serotype_parsed_results = \
        src.blastFunctions.parse_blast_results(
            args,
            serotype_output_file,
            src.genomeFunctions.get_parsing_dict('serotype')
        )

    parsed_results = None
    if serotype_parsed_results:
        parsed_results = serotype_parsed_results

    #print the requested format
    if args.tabular:
        log.info("Printing results in tabular format")

        writer = csv.writer(sys.stdout, delimiter="\t")
        writer.writerows(src.resultsToTable.results_dict_to_table(all_genomes_files
                                                                  ,all_genomes_list
                                                                  ,parsed_results))
    else:
        log.info("Printing results in JSON format")
        print(json.dumps(parsed_results, indent=4, separators=(',', ': ')))

    log.info("Done")
