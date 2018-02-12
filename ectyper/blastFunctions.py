#!/usr/bin/env python

"""
Functions for setting up, running, and parsing blast
"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import open
from builtins import str
from future import standard_library
standard_library.install_aliases()
import logging
import os

from ectyper import subprocess_util

LOG = logging.getLogger(__name__)


def create_blast_db(filelist, temp_dir):
    """http://stackoverflow.com/questions/23944657/typeerror-method-takes-1-positional-argument-but-2-were-given
    Creating a blast DB using the makeblastdb command.
    The database is created in the temporary folder of the system.

    Args:
        filelist: list of genomes that was given by the user on the command line.
        temp_dir: temporary directory to store the blastdb in.

    Returns:
        Full path to the DB
    """
    blast_db_path = os.path.join(temp_dir, 'ectyper_blastdb')

    LOG.debug("Generating the blast db at {0}".format(blast_db_path))
    cmd = [
        "makeblastdb",
        "-in", ' '.join(filelist),
        "-dbtype", "nucl",
        "-title", "ectyper_blastdb",
        "-out", blast_db_path]
    subprocess_util.run_subprocess(cmd)

    return blast_db_path


def run_blast(query_file, blast_db, args):
    """
    Execute a blastn run given the query files and blastdb

    Args:
        query_file (str): one or both of the VF / Serotype input files
        blast_db (str): validated fasta files from the user, in DB form
        args (Namespace object): parsed commadnline options from the user
        chunck_size: number of genomes in the database

    Returns:
        The blast output file
    """
    percent_identity = args.percentIdentity
    percent_length = args.percentLength

    LOG.debug('Running blast query {0} against database {1} '.format(
        query_file, blast_db))

    blast_output_file = blast_db + '.output'

    cmd = [
        "blastn",
        "-query", query_file,
        "-db", blast_db,
        "-out", blast_output_file,
        '-perc_identity', str(percent_identity),
        '-qcov_hsp_perc', str(percent_length),
        '-max_hsps', '1', # each allele only need to hit once
        # use default max_target_seqs=500
        "-outfmt",
        '6 qseqid qlen sseqid length pident sstart send sframe qcovhsp',
        "-word_size", "11"
    ]
    subprocess_util.run_subprocess(cmd)
    with open(blast_output_file, mode='rb') as fh:
        for line in fh:
            LOG.debug(line.decode('ascii'))
    return blast_output_file

def run_blast_for_identification(query_file, blast_db):
    """
    Execute a blastn run given the query files and blastdb
    with special configuration for high performance identification

    Args:
        query_file: one or both of the VF / Serotype input files
        blast_db: validated fasta files from the user, in DB form

    Returns:
        blast_output_file (str): path to the blast output file
    """

    LOG.debug('Running blast query {0} against database {1} '.format(
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
