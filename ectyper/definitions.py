#!/usr/bin/env python

"""
    Definitions for the ectyper project
"""

import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT_DIR, 'Data')
WORKPLACE_DIR = os.getcwd()

SEROTYPE_ALLELE_JSON = os.path.join(DATA_DIR, 'ectyper_alleles_db.json')
#ECOLI_MARKERS = os.path.join(DATA_DIR, 'ecoli_specific_markers.fasta')
REFSEQ_SUMMARY = os.path.join(DATA_DIR, 'assembly_summary_refseq.txt')
OSEROTYPE_GROUPS_DICT = {'1': ['O20','O137'],
                         '2': ['O28','O42'],
                         '3': ['O118','O151'],
                         '4': ['O90', 'O127'],
                         '5': ['O123', 'O186'],
                         '6': ['O46', 'O134'],
                         '7': ['O2','O50'],
                         '8': ['O107','O117'],
                         '9': ['O17','O44','O73', 'O77', 'O106'],
                         '10':['O13','O129','O135'],
                         '11':['O153','O178'],
                         '12':['O18ab', 'O18ac'],
                         '13':['O124','O164'],
                         '14':['O62','O68'],
                         '15':['O89','O101','O162'],
                         '16':['O169','O183']
                         }
