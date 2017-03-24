#!/usr/bin/env python

"""
    Definitions for the ectyper project
"""

import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
LOGGER_CONFIG = os.path.join(ROOT_DIR, 'logging.conf')
DATA_DIR = os.path.join(ROOT_DIR, 'Data/')
VF_FILE = os.path.join(DATA_DIR, 'ectyper_vfs_shortnames.ffn')
SEROTYPE_FILE = os.path.join(DATA_DIR, 'EcOH.fasta')
SEROTYPE_AND_VF_FILE = os.path.join(DATA_DIR, 'serotype_and_vfs.fasta')