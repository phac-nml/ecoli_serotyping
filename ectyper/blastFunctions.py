#!/usr/bin/env python

"""
Functions for setting up, running, and parsing blast
"""
import logging
import os

from ectyper import subprocess_util

LOG = logging.getLogger(__name__)


def create_blast_db(filelist, temp_dir):
    """
    http://stackoverflow.com/questions/23944657/typeerror-method-takes-1-positional-argument-but-2-were-given
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
