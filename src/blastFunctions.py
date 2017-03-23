#!/usr/bin/env python

"""
Functions for setting up, running, and parsing blast
"""
import collections
import subprocess
import src.genomeFunctions
import os
import tempfile
import logging

log = logging.getLogger(__name__)


def record_passes_cutoffs(blast_record, args):
    """
    For serotype prediction, need to ensure the blast hit is equal to or 
    greater than the cutoffs supplied

    :param blast_record: = {'qseqid':la[0],
                        'qlen':la[1],
                        'sseqid':la[2],
                        'length':la[3],
                        'pident':la[4]
                        }
    :param args: argparse commandline options
    :return: {True|False}
    """

    if (float(blast_record['qlen']) / float(blast_record[
                                                'length']) * 100) >= float \
                (args.percentLength) and (
                float(blast_record['pident']) >= float(args.percentIdentity)):
        return True
    else:
        return False


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
        genome_name = src.genomeFunctions.get_genome_name(
            blast_record['sseqid'])

        for blast_parser in parsing_functions:
            # we only want to store a value for a key if it doesn't exist
            # this is because blast list results from "best" to "worst" order
            # and we could be overwriting a "better" result if we do not check
            # https://www.python.org/dev/peps/pep-0448/

            parser_result = blast_parser(blast_record, args)
            results_dict[genome_name] = {**parser_result,
                                         **results_dict[genome_name]}

            # exit()
