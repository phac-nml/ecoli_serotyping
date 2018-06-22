#!/usr/bin/env python
"""
    Predictive serotyping for _E. coli_.
"""
import os
import tempfile
import datetime
import json
import logging

from ectyper import (commandLineOptions, definitions, speciesIdentification, loggingFunctions,
                     genomeFunctions, predictionFunctions, subprocess_util, __version__)

# setup the application logging
LOG = loggingFunctions.create_logger()

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
    output_directory = create_output_directory(args.output)

    #Create a file handler for log messages in the output directory
    fh = logging.FileHandler(os.path.join(output_directory, 'ectyper.log'))
    fh.setLevel(logging.DEBUG)
    LOG.addHandler(fh)

    LOG.debug(args)
    LOG.info("Starting ectyper v{}.\nOutput directory is: {}"
         .format(__version__, output_directory))

    # Initialize ectyper directories and temp files for the scope of this program
    with tempfile.TemporaryDirectory() as temp_dir:
        # Download refseq file if species identification is enabled
        if args.species:
            genomeFunctions.download_refseq()

        LOG.info("Gathering genome files")
        raw_genome_files = genomeFunctions.get_files_as_list(args.input)

        LOG.info("Identifying genome file types")
        # ['fasta']=[], ['fastq']=[]
        raw_files_dict = genomeFunctions.get_raw_files(raw_genome_files)
        LOG.debug(raw_files_dict)

        # Create files for ectyper run
        ectyper_files = create_ectyper_files(temp_dir,
                                             raw_files_dict['fastq'],
                                             output_directory
                                          )

        # Assemble any fastq files
        all_fasta_files = genomeFunctions.assembleFastq(raw_files_dict,
                                                        temp_dir,
                                                        ectyper_files['combined_fasta'],
                                                        ectyper_files['bowtie_base'])

        # Verify _E. coli_ genomes, if desired
        v_fasta_files = speciesIdentification.verify_ecoli(all_fasta_files, args.verify)

        LOG.info("Standardizing the genome headers")
        final_fasta_files = genomeFunctions.get_genome_names_from_files(v_fasta_files, temp_dir)
        LOG.info(final_fasta_files)

        # Main prediction function
        predictions_data_frame = run_prediction(final_fasta_files,
                                          args,
                                          ectyper_files['alleles_fasta'],
                                          ectyper_files['output_file'])

        # Add empty rows for genomes without a blast result
        predictions_file = predictionFunctions.add_non_predicted(
            all_fasta_files, predictions_data_frame, ectyper_files['output_file'])

        # Store most recent result in working directory
        LOG.info('\nReporting result...')
        predictionFunctions.report_result(predictions_file)


def create_output_directory(output_dir):
    """
    Create the output directory for ectyper

    :param output_dir: The user-specified output directory, if any
    :return: The output directory
    """
    # If no output directory is specified for the run, create a one based on time
    out_dir = None

    if output_dir is None:
        date_dir = ''.join([
            'ectyper_',
            str(datetime.datetime.now().date()),
            '_',
            str(datetime.datetime.now().time()).replace(':', '.')
        ])
        out_dir = os.path.join(definitions.WORKPLACE_DIR, date_dir)
    else:
        if os.path.isabs(output_dir):
            out_dir = output_dir
        else:
            out_dir = os.path.join(definitions.WORKPLACE_DIR, output_dir)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    return out_dir


def create_ectyper_files(temp_dir, fastq_list, out_dir):
    """
    Create the files needed for an ectyper run.
    This includes the fasta files and databases, and the output files.
    :param temp_dir: the temporary directory for the ectyper run
    :param fastq_list: list of all fastq files, if any
    :param out_dir: program output directory if specified
    :return: Dictionary of files for program run
    """

    # Finalize the tmp_files dictionary
    files_and_dirs = {}
    files_and_dirs['output_file'] = os.path.join(out_dir, 'output.csv')
    files_and_dirs['output_dir'] = out_dir
    files_and_dirs['alleles_fasta'] = create_alleles_fasta_file(temp_dir)
    files_and_dirs['combined_fasta'] = \
        genomeFunctions.create_combined_alleles_and_markers_file(files_and_dirs['alleles_fasta'], temp_dir)

    # Create a bowtie2 reference if there are assemblies that need to be done
    if fastq_list:
        files_and_dirs['bowtie_base'] = genomeFunctions.create_bowtie_base(temp_dir, files_and_dirs['combined_fasta'])
    else:
        files_and_dirs['bowtie_base'] = None

    LOG.info("ectyper files and directories created")
    LOG.debug(files_and_dirs)
    return files_and_dirs


def create_alleles_fasta_file(temp_dir):
    """
    Every run, re-create the fasta file of alleles to ensure a single
    source of truth for the ectyper data -- the JSON file.

    :temp_dir: temporary directory for length of program run
    :return: the filepath for alleles.fasta
    """
    output_file = os.path.join(temp_dir, 'alleles.fasta')

    with open(definitions.SEROTYPE_ALLELE_JSON, 'r') as jsonfh:
        json_data = json.load(jsonfh)

        with open(output_file, 'w') as ofh:
            for a in ["O", "H"]:
                for k in json_data[a].keys():
                    ofh.write(">" + k + "\n")
                    ofh.write(json_data[a][k]["seq"] + "\n")

    LOG.debug(output_file)
    return output_file


def run_prediction(genome_files, args, alleles_fasta, predictions_file):
    """
    Serotype prediction of all the input files, which have now been properly
    converted to fasta if required, and their headers standardized

    :param genome_files: List of genome files in fasta format
    :param args: program arguments from the commandline
    :param alleles_fasta: fasta format file of the ectyper O- and H-alleles
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
                '-query', alleles_fasta,
                '-db', blast_db,
                '-out', blast_output_file,
                '-perc_identity', str(args.percentIdentity),
                '-qcov_hsp_perc', str(args.percentLength),
                '-max_hsps', "1",
                '-outfmt', "6 qseqid qlen sseqid length pident sstart send sframe qcovhsp sseq",
                '-word_size', "11"
            ]
            subprocess_util.run_subprocess(bcline)

            LOG.info("Start serotype prediction for database #{0}".format(index + 1))
            predictions_data_frame = predictionFunctions.predict_serotype(
                blast_output_file, definitions.SEROTYPE_ALLELE_JSON, args.detailed)
        return predictions_data_frame


