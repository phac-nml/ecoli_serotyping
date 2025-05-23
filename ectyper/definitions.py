#!/usr/bin/env python

"""
    Definitions for the ectyper project
"""

import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT_DIR, 'Data')
WORKPLACE_DIR = os.getcwd()

SEROTYPE_ALLELE_JSON = os.path.join(DATA_DIR, 'ectyper_alleles_db.json')
PATHOTYPE_ALLELE_JSON = os.path.join(DATA_DIR, 'ectyper_patho_stx_toxin_typing_database.json')
SPECIES_ID_SKETCH = os.path.join(DATA_DIR, 'EnteroRef_GTDBSketch_20231003_V2.msh')
#ECOLI_MARKERS = os.path.join(DATA_DIR, 'ecoli_specific_markers.fasta')
#REFSEQ_SUMMARY = os.path.join(DATA_DIR, 'assembly_summary_refseq.txt')
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
MASH_URLS = ["https://zenodo.org/records/13969103/files/EnteroRef_GTDBSketch_20231003_V2.msh?download=1"]

HIGH_SIMILARITY_THRESHOLD_O = 0.00771 # alleles that are 99.23% apart will be reported as mixed call ~ 8 nt difference on average
MIN_O_IDENTITY_LS = 95 #low similarity group O antigen min identity threshold to pre-filter BLAST output  (identical to global threshold)
MIN_O_COVERAGE_LS = 48 #low similarity group O antigen min coverage threshold to pre-filter BLAST output (based on cross-talk study results)
PATHOTYPE_TOXIN_FIELDS = ['pathotype', 'pathotype_count', 'pathotype_genes', 'pathotype_gene_names',  'pathotype_accessions', 'pathotype_allele_id', 
                   'pathotype_pident', 'pathotype_pcov','pathotype_length_ratio', 'pathotype_rule_ids', 'pathotype_gene_counts', 'pathotype_database',
                   'stx_genes', 'stx_accessions', 'stx_allele_ids', 'stx_genes_full_name', 'stx_pidents', 'stx_pcovs', 'stx_gene_lengths', 'stx_contigs', 'stx_gene_ranges']
OUTPUT_TSV_HEADER = ['Name','Species', 'SpeciesMashRatio', 'SpeciesMashDistance','SpeciesMashTopID','O-type','H-type','Serotype','QC',
             'Evidence','GeneScores','AlleleKeys','GeneIdentities(%)',
             'GeneCoverages(%)','GeneContigNames','GeneRanges',
             'GeneLengths','DatabaseVer','Warnings','Pathotype', 'PathotypeCounts', 'PathotypeGenes', 'PathotypeGeneNames', 'PathotypeAccessions', 'PathotypeAlleleIDs', 
             'PathotypeIdentities(%)','PathotypeCoverages(%)','PathotypeGeneLengthRatios','PathotypeRuleIDs', 'PathotypeGeneCounts', 'PathoDBVer',
             'StxSubtypes','StxAccessions','StxAlleleIDs','StxAlleleNames', 'StxIdentities(%)','StxCoverages(%)','StxLengths',
             'StxContigNames','StxCoordinates']
OUTPUT_FILES_LIST = ['blastn_output_alleles.txt', 'blastn_pathotype_alleles_overall.txt', 'mash_output.txt', 
                     'stx1_allhits_annotated_df.txt', 'stx2_allhits_annotated_df.txt']