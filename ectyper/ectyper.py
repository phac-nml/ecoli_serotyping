#!/usr/bin/env python
"""
    Predictive serotyping for _E. coli_.
"""
import logging
import os
import tempfile
import datetime

from ectyper import (commandLineOptions, definitions, speciesIdentification,
                     genomeFunctions, loggingFunctions, predictionFunctions, subprocess_util)

LOG_FILE = loggingFunctions.initialize_logging()
LOG = logging.getLogger(__name__)


def run_program():
    """
    Wrapper for serotyping
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
    LOG.debug(args)
    LOG.info("Starting ectyper.\nLog file is: {}".format(LOG_FILE))

    # Initialize ectyper directories and temp files for the scope of this program
    with tempfile.TemporaryDirectory() as temp_dir:
        ectyper_files = create_ectyper_files(temp_dir, output_dir=args.output)

        # Download refseq file if species identification is enabled
        if args.species:
            genomeFunctions.download_refseq()

        LOG.info("Gathering genome files")
        raw_genome_files = genomeFunctions.get_files_as_list(args.input)

        LOG.info("Identifying genome file types")
        # ['fasta']=[], ['fastq']=[]
        raw_files_dict = genomeFunctions.get_raw_files(raw_genome_files)
        LOG.debug(raw_files_dict)

        # Assemble any fastq files
        all_fasta_files = genomeFunctions.assembleFastq(raw_files_dict, temp_dir)

        # Verify _E. coli_ genomes, if desired
        v_fasta_files = speciesIdentification.verify_ecoli(all_fasta_files, args.verify)

        LOG.info("Standardizing the genome headers")
        final_fasta_files = genomeFunctions.get_genome_names_from_files(v_fasta_files, temp_dir)
        LOG.info(final_fasta_files)

        # Main prediction function
        predictions_file = run_prediction(final_fasta_files, args,
                                          ectyper_files['output_file'])

        # Add empty rows for genomes without blast result
        predictions_file = predictionFunctions.add_non_predicted(
            raw_genome_files, predictions_file)
        LOG.info('Output saved to {0}'.format(ectyper_files['output_dir']))

        # Store most recent result in working directory
        LOG.info('\nReporting result...')
        predictionFunctions.report_result(predictions_file)


def create_ectyper_files(temp_dir, output_dir=None):
    """
    Create dictionary of files and directories used by ectyper

    Args:
        temp_dir: program scope temporary directory
        output_dir(str, optional):
            directory to store output

    Return:
        a dictionary of temporary files
        example:
            {'assemble_temp_dir': 'test/temp/assemblies',
             'fasta_temp_dir': 'test/temp/fastas',
             'output_dir': os.path.abspath('output')+'/',
             'output_file': os.path.abspath('output/output.csv')}
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

    LOG.info("ectyper files and directories created")
    LOG.debug(files_and_dirs)
    return files_and_dirs


def run_prediction(genome_files, args, predictions_file):
    """
    Serotype prediction of all the input files, which have now been properly
    converted to fasta if required, and their headers standardized

    :param genome_files: List of genome files in fasta format
    :param args: program arguments from the commandline
    :param predictions_file: the output file to store the predictions
    :return: predictions_file
    """

    # create a temp dir for blastdb
    with tempfile.TemporaryDirectory() as temp_dir:
        # Divide genome files into groups and create databases for each set
        predictions_output_file = None
        group_size = definitions.GENOME_GROUP_SIZE
        genome_groups = [
            genome_files[i:i + group_size]
            for i in range(0, len(genome_files), group_size)
        ]

        for index, g_group in enumerate(genome_groups):

            LOG.info("Creating blast database #{0} from {1}".format(index + 1, g_group))
            blast_db = os.path.join(temp_dir, "blastdb_" + str(index))
            blast_db_cmd = [
                "makeblastdb",
                "-in", ' '.join(g_group),
                "-dbtype", "nucl",
                "-title", "ectyper_blastdb",
                "-out", blast_db]
            subprocess_util.run_subprocess(blast_db_cmd)

            LOG.info("Starting blast alignment on database #{0}".format(index + 1))
            blast_output_file = blast_db + ".output"
            bcline = [
                'blastn',
                '-query', definitions.SEROTYPE_FILE,
                '-db', blast_db,
                '-out', blast_output_file,
                '-perc_identity', str(args.percentIdentity),
                '-qcov_hsp_perc', str(args.percentLength),
                '-max_hsps', "1",
                '-outfmt', "6 qseqid qlen sseqid length pident sstart send sframe qcovhsp",
                '-word_size', "11"
            ]
            subprocess_util.run_subprocess(bcline)

            LOG.info("Start serotype prediction for database #{0}".format(index + 1))
            predictions_output_file = predictionFunctions.predict_serotype(
                blast_output_file, definitions.SEROTYPE_ALLELE_JSON, predictions_file,
                args.detailed)
        return predictions_output_file


