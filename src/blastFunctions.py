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

    # We want to ensure that a match greater than 100 (due to gaps) is treated
    # as not being greater than a perfect match
    # Either direction from 100% id should be treated the same
    init_value = float(blast_record['length']) / float(blast_record['qlen']) * 100
    lid_value = float((abs(100 - init_value) * -1) % 100)

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
    completed_process = subprocess.run(["makeblastdb",
                                        "-in", ' '.join(filelist),
                                        "-dbtype", "nucl",
                                        "-title", "ectyper_blastdb",
                                        "-out", blast_db_path],
                                       check=True,
                                       universal_newlines=True,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)

    if completed_process.returncode == 0:
        log.debug("Output from makeblastdb:")
        log.debug(completed_process.stdout)
        log.debug(completed_process.stderr)
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

    completed_process = \
        subprocess.run(["blastn",
                        "-query", query_file,
                        "-db", blast_db,
                        "-out", blast_output_file,
                        "-outfmt",
                        '6 " qseqid qlen sseqid length pident sstart send sframe "',
                        "-word_size", "11"],
                       check=True,
                       universal_newlines=True,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE
                       )
    if completed_process.returncode == 0:
        log.debug("Output from blastn:")
        log.debug(completed_process.stdout)
        log.debug(completed_process.stderr)
        return blast_output_file
    else:
        log.fatal("blastn did not run successfully.\n%s",
                  completed_process.stderr)
        exit(1)


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
    result_handle = open(blast_results_file, 'r')
    results_dict = collections.defaultdict(dict)

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
        # genome name to store the parsed blast_result in
        genome_name = src.genomeFunctions.get_genome_name(
            blast_record['sseqid'])
        log.debug(genome_name)

        for blast_parser_dict in parsing_dict:
            # we only want to store a value for a key if it doesn't exist
            # this is because blast list results from "best" to "worst" order
            # and we could be overwriting a "better" result if we do not check
            # https://www.python.org/dev/peps/pep-0448/

            # the returned dict has a key specifying its type, eg. 'vf',
            # 'serotype' etc.

            blast_result_dict = blast_parser_dict['parser'](blast_record, args)

            parser_type = None
            for parser_type in blast_result_dict.keys():
                parser_results = blast_result_dict[parser_type]

            if parser_type is None:
                continue
            elif parser_type in results_dict[genome_name]:
                results_dict[genome_name][parser_type] = \
                    {**parser_results,
                     **results_dict[genome_name][parser_type]}
            else:
                results_dict[genome_name][parser_type] = parser_results

    #final prediction now that we have a dictionary of parsed results
    for blast_parser_dict in parsing_dict:
        results_dict = blast_parser_dict['predictor'](results_dict)

    return results_dict

