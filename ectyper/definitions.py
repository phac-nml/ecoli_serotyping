#!/usr/bin/env python

"""
    Definitions for the ectyper project
"""

import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT_DIR, 'Data')
WORKPLACE_DIR = os.getcwd()

GENOME_GROUP_SIZE = 50
SEROTYPE_ALLELE_JSON = os.path.join(DATA_DIR, 'ectyper_new.json')
ECOLI_MARKERS = os.path.join(DATA_DIR, 'ecoli_specific_markers.fasta')
REFSEQ_SUMMARY = os.path.join(DATA_DIR, 'assembly_summary_refseq.txt')
