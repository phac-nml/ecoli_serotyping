#!/usr/bin/env python

import src.blastFunctions
import re
import logging
import src.virulencePrediction

def test_vir_parse():
    assert parse_virulence_factors({'qseqid': 'gnd|1_O26', 'sstart': '2103817', 'pident': '94.96', 'sframe': '1', 'qlen': '643', 'send': '2103183', 'length': '635', 'sseqid': 'LT615373.1'},
                                   Namespace(csv=False, input='GCA_900096815.1_Ecoli_AG100_Sample2_M9_Assembly_genomic.fna', minimumGenomes=1, percentIdentity=90, percentLength=90, serotyper=True, virulenceFactors=True)) == \
           {'vf': {'agn43': {'description': 'antigen 43', 'blast_record': {'qlen': '3129', 'length': '3135', 'sseqid': 'LT615373.1', 'pident': '90.59', 'sframe': '1', 'sstart': '2074471', 'qseqid': 'VFG033837(gi_26249490)_(agn43)_antigen_43_[Antigen_43,_AIDA-I_type_(CVF451)]_[Escherichia_coli_CFT073]', 'send': '2077594'}}}}