#!/usr/bin/env python
"""
    Predictive serotyping for _E. coli_.
"""
import logging
import os
import sys
import tempfile
import datetime
from urllib.request import urlretrieve

from ectyper import (blastFunctions, commandLineOptions, definitions,
                     genomeFunctions, loggingFunctions, predictionFunctions,
                     speciesIdentification)

LOG_FILE = loggingFunctions.initialize_logging()
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
    (3) Query for serotype
    (4) Parse the results
    (5) Display the results
    :return: success or failure

    """
    # Initialize the program
    args = commandLineOptions.parse_command_line()
    LOG.info('Starting ectyper\nSerotype prediction with input:\n \
            {0}\n \
            Log file is: {1}'.format(args, LOG_FILE))
    LOG.debug(args)

    ## Initialize temporary directories for the scope of this program
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_files = create_tmp_files(temp_dir, output_dir=args.output)
        LOG.debug(temp_files)

        # Download refseq file if speicies identification is enabled
        if args.species:
            download_refseq()

        LOG.info("Gathering genome files")
        raw_genome_files = genomeFunctions.get_files_as_list(args.input)
        LOG.debug(raw_genome_files)

        LOG.info("Removing invalid file types")
        raw_files_dict = get_raw_files(raw_genome_files)
        LOG.debug(raw_files_dict)

        # Assembling fastq + verify ecoli genome
        LOG.info("Preparing genome files for blast alignment")
        final_fasta_files = filter_for_ecoli_files(
            raw_files_dict, temp_files, verify=args.verify, species=args.species
        )
        LOG.debug(final_fasta_files)
        if len(final_fasta_files) is 0:
            LOG.info("No valid genome files. Terminating the program.")
            exit(0)

        LOG.info("Standardizing the genome headers")
        (all_genomes_list, all_genomes_files) = \
            genomeFunctions.get_genome_names_from_files(
                final_fasta_files, temp_files['fasta_temp_dir'])
        LOG.debug(all_genomes_list)
        LOG.debug(all_genomes_files)

        # Main prediction function
        predictions_file = run_prediction(all_genomes_files, args,
                                          temp_files['output_file'])

        # Add empty rows for genomes without blast result
        predictions_file = predictionFunctions.add_non_predicted(
            all_genomes_list, predictions_file)
        LOG.info('Outputs stored in {0}'.format(temp_files['output_dir']))

        # Store most recent result in working directory
        LOG.info('\nReporting result...')
        predictionFunctions.report_result(predictions_file)

def download_refseq():
    '''
    Download refseq file with progress bar
    '''
    def download_file(url, dst):
        '''
        download file with progress bar
        '''
        urlretrieve(url, dst, reporthook)


    def reporthook(blocknum, blocksize, totalsize):
        '''
        https://stackoverflow.com/questions/15644964/python-progress-bar-and-downloads
        '''
        readsofar = blocknum * blocksize
        if totalsize > 0:
            s = "\r {:5.1%} {:{}d} / {:d}".format(
                readsofar/totalsize, readsofar,
                len(str(totalsize)),
                totalsize
            )
            sys.stderr.write(s)
            if readsofar >= totalsize: # near the end
                sys.stderr.write("\n")
        else: # total size is unknown
            sys.stderr.write("read {}\n".format(readsofar))

    if not os.path.isfile(definitions.REFSEQ_SKETCH):
        refseq_url = 'https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh'
        LOG.info("No refseq found. Downloading reference file for species identification...")
        download_file(refseq_url, definitions.REFSEQ_SKETCH)
        LOG.info("Download complete.")

def create_tmp_files(temp_dir, output_dir=None):
    """
    Return a dictionary of temporary files used by ectyper
    """

    # Get the correct files and directories
    files_and_dirs = {
        'assemble_temp_dir': os.path.join(temp_dir, 'assemblies'),
        'fasta_temp_dir': os.path.join(temp_dir, 'fastas'),
    }
    if output_dir is None:
        output_dir = ''.join([
            str(datetime.datetime.now().date()),
            '_',
            str(datetime.datetime.now().time()).replace(':', '.')
        ])
    if os.path.isabs(output_dir):
        pass
    else:
        output_dir = os.path.join(definitions.WORKPLACE_DIR, 'output', output_dir)

    output_file = os.path.join(output_dir, 'output.csv')
    if os.path.isfile(output_file):
        os.remove(output_file)
    for d in [
            output_dir, files_and_dirs['assemble_temp_dir'],
            files_and_dirs['fasta_temp_dir']
    ]:
        if not os.path.exists(d):
            os.makedirs(d)

    # Finalize the tmp_files dictionary
    files_and_dirs['output_file'] = output_file
    files_and_dirs['output_dir'] = output_dir

    LOG.info("Temporary files and directory created")
    return files_and_dirs


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
        # Divide genome files into chunks
        chunk_size = 50
        genome_chunks = [
            genome_files[i:i + chunk_size]
            for i in range(0, len(genome_files), chunk_size)
        ]
        for index, chunk in enumerate(genome_chunks):
            LOG.info("Start creating blast database #{0}".format(index + 1))
            blast_db = blastFunctions.create_blast_db(chunk, temp_dir)

            LOG.info("Start blast alignment on database #{0}".format(index + 1))
            blast_output_file = blastFunctions.run_blast(query_file, blast_db, args)
            LOG.info("Start serotype prediction for database #{0}".format(index + 1))
            predictions_file = predictionFunctions.predict_serotype(
                blast_output_file, ectyper_dict_file, predictions_file,
                args.verbose)
        return predictions_file


def get_raw_files(raw_files):
    """
    Take all the raw files, and filter not fasta / fastq

    :param raw_files:
    :return (raw_fasta_files, raw_fastq_files)
    """
    fasta_files = []
    fastq_files = []

    for file in raw_files:
        file_format = genomeFunctions.get_valid_format(file)
        if file_format == 'fasta':
            fasta_files.append(file)
        elif file_format == 'fastq':
            fastq_files.append(file)

    LOG.debug('raw fasta files: {}'.format(fasta_files))
    LOG.debug('raw fastq files: {}'.format(fastq_files))

    return({'fasta':fasta_files, 'fastq':fastq_files})


def filter_for_ecoli_files(raw_dict, temp_files, verify=False, species=False):
    """
    :param raw_dict{fasta:list_of_files, fastq:list_of_files}:
    :parapm temp_file:
    :param verify:
    :param species:
    """
    final_files = []
    for f in raw_dict.keys():
        temp_dir = temp_files['fasta_temp_dir'] if f == "fasta" else temp_files['assemble_temp_dir']

        for ffile in raw_dict[f]:
            filtered_file = filter_file_by_species(
                ffile, f, temp_dir, verify=verify, species=species)
            if filtered_file is not None and \
            genomeFunctions.get_valid_format(filtered_file) is not None:
                final_files.append(filtered_file)

    LOG.info('{} final fasta files'.format(len(final_files)))
    return final_files

def filter_file_by_species(genome_file, genome_format, temp_dir, verify=False, species=False):
    '''
    Core species recognition functionality
    :param genome_file:
    :param genome_format:
    :param temp_dir:
    :param verify:
    :param mash:
    :returns filtered_file
    '''
    combined_file = definitions.COMBINED
    filtered_file = None
    if genome_format == 'fastq':
        iden_file, pred_file = \
            genomeFunctions.assemble_reads(genome_file, combined_file, temp_dir)
        # If no alignment resut, the file is definitely not E.Coli
        if genomeFunctions.get_valid_format(iden_file) is None:
            LOG.warning(
                "{} is filtered out because no identification alignment found".format(genome_file))
            return filtered_file
        if not (verify or species) or speciesIdentification.is_ecoli_genome(
                iden_file, genome_file, mash=species):
            # final check before adding the alignment for prediction
            if genomeFunctions.get_valid_format(iden_file) != 'fasta':
                LOG.warning(
                    "{0} is filtered out because no prediction alignment found".format(genome_file))
                return filtered_file
            filtered_file = pred_file
    if genome_format == 'fasta':
        if not (verify or species) \
        or speciesIdentification.is_ecoli_genome(genome_file, mash=species):
            filtered_file = genome_file
    return filtered_file
