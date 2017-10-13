#!/usr/bin/env python

"""
Functions for setting up, running, and parsing blast
"""
import collections
import logging
import os
import tempfile

from ectyper import genomeFunctions
from ectyper import subprocess_util

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

    # We want to ensure that a match greater than 100 (due to gaps) is treated
    # as not being greater than a perfect match
    # Either direction from 100% id should be treated the same
    diff = abs(float(blast_record['qlen']) - float(blast_record['length']))
    lid_value = float(blast_record['qlen']) - diff

    if (lid_value >= float(args.percentLength)) and \
            (float(blast_record['pident']) >= float(args.percentIdentity)):
        return True
    else:
        return False


def create_blast_db(filelist):
    """http://stackoverflow.com/questions/23944657/typeerror-method-takes-1-positional-argument-but-2-were-given
    Creating a blast DB using the makeblastdb command.
    The database is created in the temporary folder of the system.

    :param filelist: genome list that was given by the user on the commandline.
    :return full path of DB
    """
    tempdir = tempfile.gettempdir()
    blast_db_path = os.path.join(tempdir, 'ectyper_blastdb')

    log.debug("Generating the blast db at %s", blast_db_path)
    cmd = [
        "makeblastdb",
        "-in", ' '.join(filelist),
        "-dbtype", "nucl",
        "-title", "ectyper_blastdb",
        "-out", blast_db_path]
    subprocess_util.run_subprocess(cmd)

    return blast_db_path


def run_blast(query_file, blast_db, args, max_genome_count=1):
    """
    Execute a blastn run given the query files and blastdb

    :param query_file: one or both of the VF / Serotype input files
    :param blast_db: validated fasta files from the user, in DB form
    :param args: parsed commandline options from the user
    :param max_genome_count: number of genome in the database [250]
    :return: the blast output file
    """
    percent_identity = args.percentIdentity
    percent_length = args.percentLength
    

    log.debug('Running blast query {0} against database {1} '.format(
        query_file, blast_db))

    blast_output_file = blast_db + '.output'

    cmd = [
        "blastn",
        "-query", query_file,
        "-db", blast_db,
        "-out", blast_output_file,
        '-perc_identity', str(percent_identity),
        '-qcov_hsp_perc', str(percent_length),
        '-max_hsps', str(max_genome_count), # at most 1 hit per genome per query
        "-outfmt",
        '6 qseqid qlen sseqid length pident sstart send sframe',
        "-word_size", "11"
    ]
    subprocess_util.run_subprocess(cmd)

    return blast_output_file

def run_blast_for_identification(query_file, blast_db):
    """
    Execute a blastn run given the query files and blastdb
    with special configuration for high performance identification

    :param query_file: one or both of the VF / Serotype input files
    :param blast_db: validated fasta files from the user, in DB form
    :return: the blast output file
    """
    
    log.debug('Running blast query {0} against database {1} '.format(
        query_file, blast_db))

    blast_output_file = blast_db + '.output'

    cmd = [
        "blastn",
        "-query", query_file,
        "-db", blast_db,
        "-out", blast_output_file,
        '-perc_identity', '90',
        '-qcov_hsp_perc', '90',
        '-max_target_seqs', '1',  # we only want to know hit/no hit
        # 10 query seq, we want at most 1 hit each
        "-outfmt",
        '6 qseqid qlen sseqid length pident sstart send sframe',
        "-word_size", "11"
    ]
    subprocess_util.run_subprocess(cmd)

    return blast_output_file


def parse_blast_results(args, blast_results_file, parsing_dict):
    """
    Given the user-defined cutoffs, return only the results that pass.
    VFs use the cutoffs directly.
    Serotype may use additional logic.

    :param args: parsed commandline options from the user
    :param blast_results_file: -outfmt 6 results of vf and/or
            serotype vs. genomes
    :param parsing_dict: functions for parsing to be applied
    :return: a dictionary of genomes and results for each
    """

    log.info("Parsing blast results in {0}".format(blast_results_file))

    result_handle = open(blast_results_file, 'r')
    results_dict = collections.defaultdict(lambda: collections.defaultdict(
        dict))

    for line in result_handle:
        clean_line = line.strip()
        la = clean_line.split()

        # We will make a dict of the results, to allow the parsing functions
        # to use the blast "name" rather than array location, which will
        # facilitate changes / additional parsers that require other info
        # later on.

        #  '6 " qseqid qlen sseqid length pident sstart send sframe "',
        blast_record = {'qseqid': la[0],
                        'qlen': la[1],
                        'sseqid': la[2],
                        'length': la[3],
                        'pident': la[4],
                        'sstart': la[5],
                        'send': la[6],
                        'sframe': la[7]
                        }

        # Initially check that the result passes the length / identity filters
        if not record_passes_cutoffs(blast_record, args):
            log.debug("The following did not pass the cutoffs:")
            log.debug(blast_record)
        else:
            # genome name to store the parsed blast_result in
            genome_name = genomeFunctions.get_genome_name(
                blast_record['sseqid'])
            log.debug(genome_name)

            # we only want to store a value for a key if it doesn't exist
            # this is because blast list results from "best" to "worst" order
            # and we could be overwriting a "better" result if we do not check
            # https://www.python.org/dev/peps/pep-0448/
            blast_result_dict = parsing_dict['parser'](blast_record,
                                                       parsing_dict['data'])

            for gene in blast_result_dict.keys():
                if gene in results_dict[genome_name][parsing_dict['type']]:
                    # test to see whether the gene is a better match
                    if new_result_is_better(
                        blast_result_dict[gene], results_dict[
                        genome_name][parsing_dict['type']][gene]
                    ):
                        results_dict[genome_name][parsing_dict['type']][gene]\
                        = \
                        blast_result_dict[gene]
                else:
                    results_dict[genome_name][parsing_dict['type']][gene] = \
                        blast_result_dict[gene]

        # final prediction now that we have a dictionary of parsed results
    result_handle.close()
    return parsing_dict['predictor'](results_dict)


def new_result_is_better(new_result, old_result):
    """
    Compare two results. If the new result is a "better" match, return TRUE.
    Otherwise return FALSE
    {'antigen': 'H11', 'blast_record': {'qseqid': '1__fliC__fliC-H11__5',
    'qlen': '1467', 'sseqid': 'JHGM01000068.1', 'length': '1459',
    'pident': '99.931', 'sstart': '29221', 'send': '27763', 'sframe': '-1'}}

    :param new_result: new blast dictionary
    :param old_result: old blast dictionary
    :return: TRUE | FALSE
    """

    new_value = float(new_result['blast_record']['length']) * float(new_result[
                                                                        'blast_record'][
                                                                        'pident'])
    old_value = float(old_result['blast_record']['length']) * float(old_result[
                                                                        'blast_record'][
                                                                        'pident'])

    if new_value > old_value:
        return 1
    else:
        return 0
