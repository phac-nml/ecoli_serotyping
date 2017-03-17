#!/usr/bin/env python

"""
    Definitions for the ectyper project
"""

import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
LOGGER_CONFIG = os.path.join(ROOT_DIR, 'logging.conf')
DATA_DIR = os.path.join(ROOT_DIR, 'Data/')
VF_FILE = os.path.join(DATA_DIR, 'repaired_ecoli_vfs_shortnames.ffn')
SEROTYPE_FILE = os.path.join(DATA_DIR, 'ECOH.fasta')