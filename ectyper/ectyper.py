#!/usr/bin/env python

"""
    Predictive genomics for _E. coli_ from the command line.
    Currently includes serotyping and VF finding.
"""

import csv
import json
import logging
import sys
import tempfile
import timeit

from ectyper import (blastFunctions, commandLineOptions, definitions,
                     genomeFunctions, loggingFunctions, resultsToTable,
                     speciesIdentification)

LOG = logging.getLogger(__name__)


def run_program():
    """
    Wrapper for both the serotyping and virulence finder
    The program needs to do the following
    (1) Filter genome files based on format
    (2) Filter genome files based on species
    (3) Map FASTQ files to FASTA seq
    (1) Get names of all genomes being tested
    (2) Create a BLAST database of those genomes
    (3) Query for serotype and/or virulence factors
    (4) Parse the results
    (5) Display the results
    :return: success or failure

    """
    start_time = timeit.default_timer()
    loggingFunctions.initialize_logging()
    args = commandLineOptions.parse_command_line()
    LOG.debug(args)
    # tempfile.tempdir is a singleton variable
    tempfile.tempdir = None
    tempdir_obj = tempfile.TemporaryDirectory()
    tempfile.tempdir = tempdir_obj.name
    # use serotype sequence from `More Serotype`
    query_file = definitions.SEROTYPE_FILE

    LOG.info("Gathering genome files")
    raw_genome_files = genomeFunctions.get_files_as_list(args.input)
    LOG.debug(raw_genome_files)

    LOG.info("Filter genome files based on format")
    raw_fasta_files = []
    raw_fastq_files = []
    for file in raw_genome_files:
        file_format = genomeFunctions.get_valid_format(file)
        if file_format == 'fasta':
            raw_fasta_files.append(file)
        elif file_format == 'fastq':
            raw_fastq_files.append(file)
    LOG.debug('raw fasta files: %s', str(raw_fasta_files))
    LOG.debug('raw fastq files: %s', str(raw_fastq_files))

    LOG.info("Filter non-ecoli genome files")
    final_fasta_files = []
    for file in raw_fasta_files:
        if speciesIdentification.is_ecoli_genome(file, args):
            final_fasta_files.append(file)
    for file in raw_fastq_files:
        iden_file, pred_file = \
            genomeFunctions.assemble_reads(file, definitions.COMBINED)
        if speciesIdentification.is_ecoli_genome(iden_file, args, file):
            final_fasta_files.append(pred_file)
    
    LOG.info('Final fasta files: %s', str(final_fasta_files))

    if final_fasta_files == []:
        LOG.info("No valid genome file. Terminating the program.")
        tempdir_obj.cleanup()
        exit(1)

    LOG.info("Gathering genome names from files")
    (all_genomes_list, all_genomes_files) = \
        genomeFunctions.get_genome_names_from_files(final_fasta_files)
    LOG.debug(all_genomes_list)
    LOG.debug(all_genomes_files)

    num_of_genome = len(all_genomes_files)

    LOG.info("Creating blast database")
    blast_db = blastFunctions.create_blast_db(all_genomes_files)

    serotype_output_file = \
        blastFunctions.run_blast(query_file, blast_db, args, num_of_genome)
    parsed_results = \
        blastFunctions.parse_blast_results(
            args,
            serotype_output_file,
            genomeFunctions.get_parsing_dict('serotype')
        )

    if args.tabular:
        LOG.info("Printing results in tabular format")

        writer = csv.writer(sys.stdout, delimiter="\t")
        writer.writerows(
            resultsToTable.results_dict_to_table(
                all_genomes_files,
                all_genomes_list,
                parsed_results))
    else:
        LOG.info("Printing results in JSON format")
        LOG.info(json.dumps(parsed_results, indent=4, separators=(',', ': ')))

    elapsed_time = timeit.default_timer() - start_time
    LOG.info("Program completed successfully in %0.3f sec.", elapsed_time)
    final_output = {}
    tempdir_obj.cleanup()

    # simplify output
    # [{'genome_name':xxx,'predicted O':xxx,'predicted H':xxx}]
    output = []
    for key, value in parsed_results.items():
        output_entry = {
            'genome name': key,
            'predicted O': None,
            'predicted H': None
        }
        try:
            output_entry['predicted O'] = value['serotype']['otype']
        except KeyError:
            pass
        try:
            output_entry['predicted H'] = value['serotype']['htype']
        except KeyError:
            pass
        output.append(output_entry)
    LOG.info('\nSummary:\n%s',json.dumps(output, indent=4, separators=(',', ': ')))
    return output