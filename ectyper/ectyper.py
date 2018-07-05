#!/usr/bin/env python
"""
    Predictive serotyping for _E. coli_.
"""
import os
import tempfile
import datetime
import json
import logging
from multiprocessing import Pool

from ectyper import (commandLineOptions, definitions, speciesIdentification, loggingFunctions,
                     genomeFunctions, predictionFunctions, subprocess_util, __version__)

# setup the application logging
LOG = loggingFunctions.create_logger()

def run_program():
    """
    Main function for E. coli serotyping.
    Creates all required files and controls function execution.
    :return: success or failure
    """

    args = commandLineOptions.parse_command_line()
    output_directory = create_output_directory(args.output)

    #Create a file handler for log messages in the output directory
    fh = logging.FileHandler(os.path.join(output_directory, 'ectyper.log'))
    fh.setLevel(logging.DEBUG)
    LOG.addHandler(fh)

    LOG.debug(args)
    LOG.info("Starting ectyper v{}\nOutput directory is: {}"
         .format(__version__, output_directory))

    # Initialize ectyper directory for the scope of this program
    with tempfile.TemporaryDirectory() as temp_dir:
        LOG.info("Gathering genome files")
        raw_genome_files = genomeFunctions.get_files_as_list(args.input)

        LOG.info("Identifying genome file types")
        # 'fasta'=[], 'fastq'=[], 'other'=[]
        raw_files_dict = genomeFunctions.get_raw_files(raw_genome_files)

        alleles_fasta = create_alleles_fasta_file(temp_dir)
        combined_fasta = genomeFunctions.create_combined_alleles_and_markers_file(alleles_fasta, temp_dir)
        bowtie_base = genomeFunctions.create_bowtie_base(temp_dir, combined_fasta) if raw_files_dict['fastq'] else None

        # Assemble any fastq files, get final fasta list
        all_fasta_files = genomeFunctions.assembleFastq(raw_files_dict,
                                                        temp_dir,
                                                        combined_fasta,
                                                        bowtie_base)

        # Verify we have at least one fasta file. Optionally species ID.
        # Get a tuple of ecoli and other genomes
        (ecoli_genomes, other_genomes_dict) = speciesIdentification.verify_ecoli(all_fasta_files,
                                                                                 raw_files_dict['other'],
                                                                                 args)

        LOG.info("Standardizing the genome headers")
        final_fasta_files = genomeFunctions.get_genome_names_from_files(ecoli_genomes, temp_dir)

        # Main prediction function
        predictions_dict = run_prediction(final_fasta_files,
                                          args,
                                          alleles_fasta)

        # Add empty rows for genomes without a blast result
        final_predictions = predictionFunctions.add_non_predicted(
            raw_genome_files, predictions_dict, other_genomes_dict)

        # Store most recent result in working directory
        LOG.info("Reporting results:\n")

        predictionFunctions.report_result(final_predictions, os.path.join(output_directory, 'output.csv'))
        LOG.info("\nECTyper has finished successfully.")


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


def run_prediction(genome_files, args, alleles_fasta):
    """
    Serotype prediction of all the input files, which have now been properly
    converted to fasta if required, and their headers standardized

    :param genome_files: List of genome files in fasta format
    :param args: program arguments from the commandline
    :param alleles_fasta: fasta format file of the ectyper O- and H-alleles
    :param predictions_file: the output file to store the predictions
    :return: predictions_file
    """

    predictions_dict = {}
    # create a temp dir for blastdb
    with tempfile.TemporaryDirectory() as temp_dir:
        # Divide genome files into groups and create databases for each set
        per_core = int(len(genome_files) / args.cores)
        group_size = 50 if per_core > 50 else per_core

        genome_groups = [
            genome_files[i:i + group_size]
            for i in range(0, len(genome_files), group_size)
        ]
        LOG.info("Creating blast databases")
        for index, g_group in enumerate(genome_groups):

            LOG.debug("Creating blast database #{0} from {1}".format(index + 1, g_group))
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

            LOG.info("Starting serotype prediction for database #{0}".format(index + 1))
            db_prediction_dict = predictionFunctions.predict_serotype(
                blast_output_file, definitions.SEROTYPE_ALLELE_JSON)

            # merge the per-database predictions with the final predictions dict
            predictions_dict = {**db_prediction_dict, **predictions_dict}
        return predictions_dict


