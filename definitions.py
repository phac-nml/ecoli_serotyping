#!/usr/bin/env python

"""
    Definitions for the ectyper project
"""

import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
LOGGER_CONFIG = os.path.join(ROOT_DIR, 'logging.conf')
DATA_DIR = os.path.join(ROOT_DIR, 'Data')
# original
SEROTYPE_FILE = os.path.join(DATA_DIR, 'EcOH.fasta')
SEROTYPE_ALLELE_JSON = os.path.join(DATA_DIR, 'allele_types.json')
# from MoreSerotype
SEROTYPE_FILE = os.path.join(DATA_DIR, 'serotype_dict.fasta')
SEROTYPE_ALLELE_JSON = os.path.join(DATA_DIR, 'allele_serotype.json')
