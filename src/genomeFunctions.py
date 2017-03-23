#!/usr/bin/env python

import os
import logging.config
import re
import Bio.SeqIO
import tempfile
import subprocess
import src.serotypePrediction
import src.virulencePrediction
import collections

log = logging.getLogger(__name__)


def get_files_as_list(file_or_directory):
    """
    Creates a list of files from either the given file, or all files within the
    directory specified (where each file name is its absolute path).

    :param file_or_directory: file or directory name given on commandline
    :return files_list: List of all the files found.

    """

    files_list = []

    if os.path.isdir(file_or_directory):
        log.info("Gathering genomes from directory " + file_or_directory)

        # Create a list containing the file names
        for root, dirs, files in os.walk(file_or_directory):
            for filename in files:
                files_list.append(os.path.join(root, filename))

    else:
        log.info("Using genomes in file " + file_or_directory)
        files_list.append(os.path.abspath(file_or_directory))

    # check that all are valid fasta files
    # if not, exclude with warning
    validated_files_list = validate_fasta_files(files_list)

    return sorted(validated_files_list)


def validate_fasta_files(files):
    """
    Check using Bio.SeqIO if files are valid fasta format.

    :param files: full path of all files
    :return: a list of all files that pass
    """

    validated_files = []
    for file in files:
        for _ in Bio.SeqIO.parse(file, "fasta"):
            log.debug("%s is a valid fasta file", file)
            validated_files.append(file)

            break

    return validated_files


def get_genome_names_from_files(files):
    """
    For each file:
    Takes the first header from a fasta file and sends it to the get_genome_name
    function. Returns the name of the genome. If the name of the file is to be
    used as the genome name, creates a temporary file using >lcl|filename as the
    first part of the header.

    :param files: The list of files to get the genome names for
    :return: ([genome names], [file names])
    """

    list_of_genomes = []
    list_of_files = []
    for file in files:
        header = get_fasta_header_from_file(file)
        genome_name = get_genome_name(header)

        # if the header and genome_name are the same, we need to use the filename
        # as the genome name. This means we also need to create a new file adding
        # the filename to each of the headers, so that downstream applications
        # (eg. BLAST) can be used with the filename as genome name.

        if header == genome_name:
            # get only the name of the file for use in the fasta header
            file_path_name = os.path.splitext(os.path.basename(file))
            n_name = file_path_name[0]

            # create a new file for the updated fasta headers
            new_file_tuple = tempfile.mkstemp()
            new_file = new_file_tuple[1]

            # add the new name to the list of files and genomes
            list_of_files.append(new_file)
            list_of_genomes.append(n_name)

            with open(new_file, "w") as output_fh:
                for record in Bio.SeqIO.parse(file, "fasta"):
                    record.description = ">lcl|" + n_name + "|" + record.description
                    log.debug(record.description)
                    Bio.SeqIO.write(record, output_fh, "fasta")
        else:
            list_of_files.append(file)
            list_of_genomes.append(genome_name)

    return list_of_genomes, list_of_files


def get_genome_name(header):
    """
    Getting the name of the genome by hierarchy. This requires reading the first
    fasta header from the file. It also assumes a single genome per file.

    :param header: The header containing the record.
    :return genomeName: Name of the genome contained in the header.
    """

    re_patterns = (
        # Look for lcl followed by the possible genome name
        re.compile('lcl\|([\w-]*)'),

        # Look for a possible genome name at the beginning of the record ID
        re.compile('^(\w{8}\.\d)'),

        # Look for ref, gb, emb or dbj followed by the possible genome name
        re.compile('(ref\|\w{2}_\w{6}|gb\|\w{8}|emb\|\w{8}|dbj\|\w{8})'),

        # Look for gi followed by the possible genome name
        re.compile('(gi\|\d{8})')
    )

    # if nothing matches, use the full header as genome_name
    genome_name = header
    for rep in re_patterns:
        m = rep.search(header)

        if m:
            genome_name = m.group(1)
            break

    return genome_name


def get_fasta_header_from_file(filename):
    """
    Gets the first fasta sequence from the file, and returns the fasta header.
    The files should have already been validated as fasta format.

    :param filename: the absolute path of the fasta file
    :return: header
    """

    for record in Bio.SeqIO.parse(filename, "fasta"):
        return record.description


def create_blast_db(filelist):
    """
    Creating a blast DB using the makeblastdb command.
    The database is created in the temporary folder of the system.

    :param filelist: genome list that was given by the user on the commandline.
    :return full path of DB
    """

    tempdir = tempfile.mkdtemp()
    blast_db_path = os.path.join(tempdir, 'ectyper_blastdb')

    log.debug("Generating the blast db at %s", blast_db_path)
    completed_process = subprocess.run(["makeblastdb",
                                        "-in", ' '.join(filelist),
                                        "-dbtype", "nucl",
                                        "-title", "ectyper_blastdb",
                                        "-out", blast_db_path],
                                       check=True)

    if completed_process.returncode == 0:
        return blast_db_path
    else:
        log.fatal("makeblastdb was unable to successfully run.\n%s",
                  completed_process.stderr)


def run_blast(query_file, blast_db):
    """
    Execute a blastn run given the query files and blastdb

    :param query_file: one or both of the VF / Serotype input files
    :param blast_db: validated fasta files from the user, in DB form
    :return: the blast output file
    """

    blast_output_file = blast_db + '.output'

    completed_process = subprocess.run(["blastn",
                                        "-query", query_file,
                                        "-db", blast_db,
                                        "-out", blast_output_file,
                                        "-outfmt",
                                        '6 " qseqid qlen sseqid length pident "',
                                        "-word_size", "11"])
    if completed_process.returncode == 0:
        return blast_output_file
    else:
        log.fatal("blastn did not run successfully.\n%s",
                  completed_process.stderr)
        exit(1)


def get_list_of_parsing_functions(args):
    """
    Given the parsed arguments from argparser, return a list of functions.
    :param args:
    :return: [function]
    """

    list_of_functions = []
    if args.serotyper:
        list_of_functions.append(src.serotypePrediction.predict_serotype)

    if args.virulenceFactors:
        list_of_functions.append(
            src.virulencePrediction.predict_virulence_factors)
    return list_of_functions


def parse_blast_results(args, blast_results_file, parsing_functions):
    """
    Given the user-defined cutoffs, return only the results that pass.
    VFs use the cutoffs directly.
    Serotype may use additional logic.

    :param args: parsed commandline options from the user
    :param blast_results_file: -outfmt 6 results of vf and/or serotype vs. genomes
    :param parsing_functions: functions for parsing to be applied
    :return: a dictionary of genomes and results for each
    """
    result_handle = open(blast_results_file, 'r')
    results_dict = collections.defaultdict(dict)

    for line in result_handle:
        clean_line = line.strip()
        la = clean_line.split()

        # We will make a dict of the results, to allow the parsing functions
        # to use the blast "name" rather than array location, which will
        # facilitate changes / additional parsers that require other info
        # later on.

        blast_record = {'qseqid': la[0],
                        'qlen': la[1],
                        'sseqid': la[2],
                        'length': la[3],
                        'pident': la[4]
                        }

        # genome name to store the parsed blast_result in

        genome_name = get_genome_name(blast_record['sseqid'])

        for blast_parser in parsing_functions:
            # we only want to store a value for a key if it doesn't exist
            # this is because blast list results from "best" to "worst" order
            # and we could be overwriting a "better" result if we do not check
            # https://www.python.org/dev/peps/pep-0448/

            results_dict[genome_name] = {
            **blast_parser(blast_record, args),
            **results_dict[genome_name]}

            #log.debug(genome_name)
            # exit()
