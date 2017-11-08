#!/usr/bin/env python

"""
    Predictive genomics for _E. coli_ from the command line.
    Currently includes serotyping and VF finding.
"""
import logging
import os
import tempfile
import timeit
import datetime
import sys

from ectyper import (blastFunctions, commandLineOptions, definitions,
                     genomeFunctions, loggingFunctions, predictionFunctions,
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
    # Initialize the program
    ## Initialize logging and timer
    start_time = timeit.default_timer()
    curr_time = timeit.default_timer()
    loggingFunctions.initialize_logging()
    ## Initialize temporary directories for the scope of this program
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create temporary folders
        assemble_temp_dir = os.path.join(temp_dir, 'assembles')
        os.mkdir(assemble_temp_dir)
        fasta_temp_dir = os.path.join(temp_dir, 'fastas')
        os.mkdir(fasta_temp_dir)
        ## Parse arguments
        args = commandLineOptions.parse_command_line(sys.argv[1:])
        LOG.info('\nStarting ectyper')
        LOG.debug(args)
        ## Get constants from definitions
        workplace_dir = definitions.WORKPLACE_DIR
        query_file = definitions.SEROTYPE_FILE
        ectyper_dict_file = definitions.SEROTYPE_ALLELE_JSON
        # Create output directory
        output_file = os.path.join(
            workplace_dir,
            'output',
            str(datetime.datetime.now().date()) + '_' + str(datetime.datetime.now().time()).replace(':', '.'),
            'output.csv')
        output_dir = os.path.split(output_file)[0]
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            LOG.info('output directory is created')

        # Collect genome files
        LOG.info("Gathering genome files")
        raw_genome_files = genomeFunctions.get_files_as_list(args.input)
        LOG.debug(raw_genome_files)

        curr_time = timeit.default_timer()
        # Filter invalid file formats
        LOG.info("Start filtering based on file format")
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
        LOG.info("Finished filtering based on file format in %.2f seconds",
                timeit.default_timer()-curr_time)
        curr_time = timeit.default_timer()

        # Filter invalid species
        LOG.info("Start filtering non-ecoli genome files")
        final_fasta_files = []
        for file in raw_fasta_files:
            filtered_file = filter_file_by_species(file, 'fasta', assemble_temp_dir)
            if filtered_file:
                final_fasta_files.append(filtered_file)
        for file in raw_fastq_files:
            filtered_file = filter_file_by_species(file, 'fastq', assemble_temp_dir)
            if filtered_file:
                final_fasta_files.append(filtered_file)
        LOG.info("Finished filtering non-ecoli genome files in %.2f seconds",
                timeit.default_timer()-curr_time)
        curr_time = timeit.default_timer()

        LOG.info('%d final fasta files', len(final_fasta_files))
        if final_fasta_files == []:
            LOG.info("No valid genome file. Terminating the program.")
            exit(1)
        # Convert genome headers
        LOG.info("Start standardize genome headers")
        (all_genomes_list, all_genomes_files) = \
            genomeFunctions.get_genome_names_from_files(final_fasta_files, fasta_temp_dir)
        LOG.debug(all_genomes_list)
        LOG.debug(all_genomes_files)
        LOG.info("Finished standardize genome headers in %.2f seconds",
                timeit.default_timer()-curr_time)
        curr_time = timeit.default_timer()

        # Main prediction function
        predictions_file = run_prediction(all_genomes_files, args, output_file)
        # Add empty rows for genomes without blast result
        predictions_file = predictionFunctions.add_non_predicted(all_genomes_list, predictions_file)

        LOG.info("Ectyper completed successfully in %0.3f sec.", timeit.default_timer() - start_time)
        LOG.info('\nReporting result...')
        predictionFunctions.report_result(predictions_file)

def run_prediction(genome_files, args, predictions_file):
    '''
    Core prediction functionality
    :param genome_files:
    :param args:
    :param predictions_file:
    :returns predictions_file
    '''
    query_file = definitions.SEROTYPE_FILE
    ectyper_dict_file = definitions.SEROTYPE_ALLELE_JSON
    # create a temp dir for blastdb
    with tempfile.TemporaryDirectory() as temp_dir:
        # Divide genome files into chunks of size 100
        chunk_size = 1000
        genome_chunks = [genome_files[i:i + chunk_size]
                        for i in range(0, len(genome_files), chunk_size)]
        for index, chunk in enumerate(genome_chunks):
            curr_time = timeit.default_timer()
            LOG.info("Start creating blast database #%d", index+1)
            blast_db = blastFunctions.create_blast_db(chunk, temp_dir)
            LOG.info("Finished creating blast database #%d in %.2f seconds",
                    index+1, timeit.default_timer()-curr_time)
            curr_time = timeit.default_timer()

            LOG.info("Start blast alignment on database #%d", index+1)
            blast_output_file = blastFunctions.run_blast(
                query_file, blast_db, args, len(chunk))
            LOG.info("Finished blast alignment on database #%d in %.2f seconds",
                    index+1, timeit.default_timer()-curr_time)
            curr_time = timeit.default_timer()
            LOG.info("Start serotype prediction for database #%d", index+1)
            predictions_file = predictionFunctions.predict_serotype(
                blast_output_file, ectyper_dict_file, predictions_file, args.verbose)
            LOG.info("Finished serotype prediction for database #%d in %.2f seconds",
                    index+1, timeit.default_timer()-curr_time)
            curr_time = timeit.default_timer()
        return predictions_file

def filter_file_by_species(genome_file, genome_format, temp_dir, mash=False):
    '''
    Core species recognition functionality
    :param genome_file:
    :param genome_format:
    :param temp_dir:
    :param mash:
    :returns filtered_file
    '''
    combined_file = definitions.COMBINED
    filtered_file = None
    if genome_format is 'fastq':
        iden_file, pred_file = \
            genomeFunctions.assemble_reads(genome_file, combined_file, temp_dir)
        # If no alignment resut, the file is definitely not E.Coli
        if genomeFunctions.get_valid_format(iden_file) is None:
            LOG.warning("%s is filtered out because no identification alignment found", genome_file)
            return filtered_file
        if speciesIdentification.is_ecoli_genome(iden_file, genome_file, mash=mash):
            # final check before adding the alignment for prediction
            if genomeFunctions.get_valid_format(iden_file) != 'fasta':
                LOG.warning("%s is filtered out because no prediction alignment found", genome_file)
                return filtered_file
            filtered_file = pred_file
    if genome_format is 'fasta':
        if speciesIdentification.is_ecoli_genome(genome_file, mash=mash):
            filtered_file = genome_file
    return filtered_file