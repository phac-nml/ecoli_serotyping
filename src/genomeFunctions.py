#!/usr/bin/env python

import os
import logging.config
import re
import Bio.SeqIO


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

    return sorted(files_list)


def get_genome_name(filename):
    """
    Getting the name of the genome by hierarchy. This requires reading the first
    fasta header from the file. It also assumes a single genome per file.

    :param filename: Name of the file containing the record.
    :return genomeName: Name of the genome contained in the file (or sequence).
    """

    record_id = get_fasta_header_from_file(filename)

    # Look for lcl followed by the possible genome name
    if re.search('lcl\|([\w-]*)', record_id):
        match = re.search('lcl\|([\w-]*)', record_id)
        match = str(match.group())
        genome_name = match.split('|')[1]

    # Look for a possible genome name at the beginning of the record ID
    elif re.search('(^[a-zA-Z][a-zA-Z]\w{6}\.\d)', record_id):
        match = re.search('(\w{8}\.\d)', record_id)
        genome_name = str(match.group())

    # Look for ref, gb, emb or dbj followed by the possible genome name
    elif re.search('(ref\|\w{2}_\w{6}|gb\|\w{8}|emb\|\w{8}|dbj\|\w{8})',
                   record_id):
        match = re.search('(ref\|\w{2}_\w{6}|gb\|\w{8}|emb\|\w{8}|dbj\|\w{8})',
                          record_id)
        match = str(match.group())
        genome_name = match.split('|')[1]

    # Look for gi followed by the possible genome name
    elif re.search('gi\|\d{8}', record_id):
        match = re.search('gi\|\d{8}', record_id)
        match = str(match.group())
        genome_name = match.split('|')[1]

    # Assign the file name as genome name
    else:
        genome_name = filename

    return genome_name


def get_fasta_header_from_file(filename):
    """
    Gets the first fasta sequence from the file, and returns the fasta header.
    Returns FALSE if the file is not a valid fasta sequence.

    :param filename: the absolute path of the fasta file
    :return: header || FALSE
    """

    for record in Bio.SeqIO.parse(filename, "fasta"):
        match = re.search('(^[a-zA-Z]+)', str(record.seq))

        if match:
            log.debug("%s is a fasta file", filename)
            return record.description
        else:
            log.debug("%s is not a fasta file", filename)
            return False