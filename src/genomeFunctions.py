#!/usr/bin/env python

import os
import logging.config
import re


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
    Getting the name of the genome by hierarchy
    :param filename: Name of the file containing the record.
    :return genomeName: Name of the genome contained in the file (or sequence).
    """

    # Look for lcl followed by the possible genome name
    if re.search('lcl\|([\w-]*)', recordID):
        match = re.search('lcl\|([\w-]*)', recordID)
        match = str(match.group())
        genome_name = match.split('|')[1]

    # Look for a possible genome name at the beginning of the record ID
    elif re.search('(^[a-zA-Z][a-zA-Z]\w{6}\.\d)', recordID):
        match = re.search('(\w{8}\.\d)', recordID)
        genome_name = str(match.group())

    # Look for ref, gb, emb or dbj followed by the possible genome name
    elif re.search('(ref\|\w{2}_\w{6}|gb\|\w{8}|emb\|\w{8}|dbj\|\w{8})',
                   recordID):
        match = re.search('(ref\|\w{2}_\w{6}|gb\|\w{8}|emb\|\w{8}|dbj\|\w{8})',
                          recordID)
        match = str(match.group())
        genome_name = match.split('|')[1]

    # Look for gi followed by the possible genome name
    elif re.search('gi\|\d{8}', recordID):
        match = re.search('gi\|\d{8}', recordID)
        match = str(match.group())
        genome_name = match.split('|')[1]
    # Assign the file name as genome name
    else:
        genome_name = filename

    if recordID not in FILENAMES:
        FILENAMES[recordID] = filename
    if filename not in GENOMENAMES:
        GENOMENAMES[filename] = recordID

    return genome_name
