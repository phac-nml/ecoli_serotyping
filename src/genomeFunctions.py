#!/usr/bin/env python

import logging.config
import re
import tempfile

import Bio.SeqIO
import os

import src.serotypePrediction
import src.virulencePrediction

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
                    output_fh.write(">lcl|" + n_name + "|" + record.description + "\n")
                    output_fh.write(str(record.seq) + "\n")
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
        re.compile('(lcl\|[\w\-\.]+)'),

        # Look for a possible genome name at the beginning of the record ID
        re.compile('^(\w{8}\.\d)'),

        # Look for ref, gb, emb or dbj followed by the possible genome name
        re.compile('(ref\|\w{2}_\w{6}|gb\|\w{8}|emb\|\w{8}|dbj\|\w{8})'),

        # Look for gi followed by the possible genome name
        re.compile('(gi\|\d{8})'),

        # Look for name followed by space, then description
        re.compile('^([\w\-\.]+)\s+[\w\-\.]+')
    )

    # if nothing matches, use the full header as genome_name
    genome_name = header
    for rep in re_patterns:
        m = rep.search(header)

        if m:
            genome_name = m.group(1)
            break

    return str(genome_name)


def get_fasta_header_from_file(filename):
    """
    Gets the first fasta sequence from the file, and returns the fasta header.
    The files should have already been validated as fasta format.

    :param filename: the absolute path of the fasta file
    :return: header
    """

    for record in Bio.SeqIO.parse(filename, "fasta"):
        return record.description


def get_list_of_parsing_dict(args):
    """
    Given the parsed arguments from argparser, return a dictionary of functions.
    :param args: Parsed commandline args
    :return: {name: {parser: function, predictor: function}}
    """

    list_of_dict = []
    if args.serotyper:
        list_of_dict.append({'parser':src.serotypePrediction.parse_serotype,
                            'predictor':src.serotypePrediction.predict_serotype
                             })

    if args.virulenceFactors:
        list_of_dict.append({
            'parser':src.virulencePrediction.parse_virulence_factors,
            'predictor':src.virulencePrediction.predict_virulence_factors})
    return list_of_dict

