#!/usr/bin/env python

"""
    Definitions for the ectyper project
"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

from future import standard_library
standard_library.install_aliases()
import os
import sys

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT_DIR, 'Data')
# Python3 vs Python2 difference.
try:
    # Python3
    WORKPLACE_DIR = os.getcwdu()
except:
    # Python2
    WORKPLACE_DIR = os.getcwd()

SEROTYPE_FILE = os.path.join(DATA_DIR, 'ectyper_data.fasta')
SEROTYPE_ALLELE_JSON = os.path.join(DATA_DIR, 'ectyper_dict.json')
COMBINED = os.path.join(DATA_DIR, 'combined.fasta')

ECOLI_MARKERS = os.path.join(DATA_DIR, 'ecoli_specific_markers.fasta')
SAMTOOLS = 'samtools'
REFSEQ_SUMMARY = os.path.join(DATA_DIR, 'assembly_summary_refseq.txt')
REFSEQ_SKETCH = os.path.join(DATA_DIR, 'refseq.genomes.k21s1000.msh')

try:
    # Python3
    from tempfile import TemporaryDirectory, NamedTemporaryFile
except:
    # Python2
    from backports.tempfile import TemporaryDirectory, NamedTemporaryFile
# Aliases
TEMPDIR = TemporaryDirectory
NAMEDTEMPFILE = NamedTemporaryFile
