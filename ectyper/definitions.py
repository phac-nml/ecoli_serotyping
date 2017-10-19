#!/usr/bin/env python

"""
    Definitions for the ectyper project
"""

import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
LOGGER_CONFIG = os.path.join(ROOT_DIR, 'logging.conf')
DATA_DIR = os.path.join(ROOT_DIR, 'Data')

SEROTYPE_FILE = os.path.join(DATA_DIR, 'ectyper_data.fasta')
SEROTYPE_ALLELE_JSON = os.path.join(DATA_DIR, 'ectyper_dict.json')
COMBINED = os.path.join(DATA_DIR, 'combined.fasta')

LEGACY_SEROTYPE_FILE = os.path.join(DATA_DIR, 'legacy_ectyper_data.fasta')
LEGACY_SEROTYPE_ALLELE_JSON = os.path.join(DATA_DIR, 'legacy_ectyper_dict.json')
LEGACY_COMBINED = os.path.join(DATA_DIR, 'legacy_combined.fasta')

REFERENCE_INDEX = os.path.join(DATA_DIR, 'bowtie_index/serotype_dict')
ECOLI_MARKERS = os.path.join(DATA_DIR, 'ecoli_specific_markers.fasta')
SAMTOOLS = 'samtools'
REFSEQ_SUMMARY = os.path.join(DATA_DIR, 'assembly_summary_refseq.txt')
REFSEQ_SKETCH = os.path.join(DATA_DIR, 'refseq.genomes.k21s1000.msh')

OUTPUT_DIR = os.path.join(ROOT_DIR, 'output')
