#!/usr/bin/env python

"""
    Predictive genomics for _E. coli_ from the command line.
    Currently includes serotyping and VF finding.
"""

import csv
import json
import logging
import sys

import definitions
import src.blastFunctions
import src.commandLineOptions
import src.genomeFunctions
import src.loggingFunctions
import src.resultsToTable

LOG = logging.getLogger(__name__)


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
    src.loggingFunctions.initialize_logging()
    args = src.commandLineOptions.parse_command_line()
    LOG.debug(args)

    # use serotype sequence from `More Serotype`
    query_file = definitions.SEROTYPE_FILE

    LOG.info("Gathering genome files")
    raw_genome_files = src.genomeFunctions.get_files_as_list(args.input)
    LOG.debug(raw_genome_files)

    LOG.info("Filter genome files based on format")
    raw_fasta_files = []
    raw_fastq_files = []
    for file in raw_genome_files:
        file_format = src.genomeFunctions.get_valid_format(file)
        if file_format == 'fasta':
            raw_fasta_files.append(file)
        elif file_format == 'fastq':
            raw_fastq_files.append(file)
    LOG.debug('raw fasta files: %s', str(raw_fasta_files))
    LOG.debug('raw fastq files: %s', str(raw_fastq_files))

    LOG.info("Filter non-ecoli genome files")
    final_fasta_files = []
    for file in raw_fasta_files:
        if src.genomeFunctions.is_ecoli_genome(file, args):
            final_fasta_files.append(file)
    final_fastq_files = []
    for file in raw_fastq_files:
        if src.genomeFunctions.is_ecoli_genome(file, args, True):
            final_fastq_files.append(file)

    LOG.debug('Final fasta files: %s', str(final_fasta_files))
    LOG.debug('Final fastq files: %s', str(final_fastq_files))

    if final_fastq_files != []:
        LOG.info("Assemble reads files by mapping to query file")
        for file in final_fastq_files:
            assemble = \
                src.genomeFunctions.assemble_reads(file, query_file)
            final_fasta_files.append(assemble)
    
    LOG.info('Final fasta files: %s', str(final_fasta_files))

    if final_fasta_files == []:
        LOG.info("No valid genome file. Terminating the program.")
        exit(1)

    LOG.info("Gathering genome names from files")
    (all_genomes_list, all_genomes_files) = \
        src.genomeFunctions.get_genome_names_from_files(final_fasta_files)
    LOG.debug(all_genomes_list)
    LOG.debug(all_genomes_files)

    LOG.info("Creating blast database")
    blast_db = src.blastFunctions.create_blast_db(all_genomes_files)

    serotype_output_file = \
        src.blastFunctions.run_blast(query_file, blast_db, args.percentIdentity)
    parsed_results = \
        src.blastFunctions.parse_blast_results(
            args,
            serotype_output_file,
            src.genomeFunctions.get_parsing_dict('serotype')
        )

    if args.tabular:
        LOG.info("Printing results in tabular format")

        writer = csv.writer(sys.stdout, delimiter="\t")
        writer.writerows(
            src.resultsToTable.results_dict_to_table(
                all_genomes_files,
                all_genomes_list,
                parsed_results))
    else:
        LOG.info("Printing results in JSON format")
        print(json.dumps(parsed_results, indent=4, separators=(',', ': ')))

    LOG.info("Done")
