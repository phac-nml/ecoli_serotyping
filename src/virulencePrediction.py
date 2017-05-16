#!/usr/bin/env python

import src.blastFunctions
import re
import logging

log = logging.getLogger(__name__)

"""
    Virulence prediction for E. coli
"""


def parse_virulence_factors(blast_record):
    """
    Entry point for virulence factor prediction
    :return: dict(
        otype = O
        htype = H
    )
    """

    # >VFG043647(gi_290434)_(clpH)_fimbrial_protein_[CS31A_capsule-like_antigen_(AI036)]_[Escherichia_coli]
    # >SPG000138_(ccdb)_Plasmid_KIL19_(from_E.coli)_cytotoxic_protein_(ccdB)_gene,_complete_cds._[L27082_1-381]


    vf_results = {}

    vf_re_patterns = (re.compile('^[^_]+_[^_]+_\((\w+)\)_([^[]+)'),
                      re.compile('^[^_]+_\((\w+)\)_([^[]+)')
                     )

    log.debug("Searching for VF in " + str(blast_record['qseqid']))
    for rep in vf_re_patterns:
        m = rep.search(blast_record['qseqid'])

        if m:
            vf_results['vf']= \
                {m.group(1): {
                    'description':m.group(2).replace('_', ' ').strip(),
                    'blast_record':blast_record}
                }
            log.debug(vf_results)
            break

    return vf_results


def predict_virulence_factors(results_dict):
    """
    
    :return: 
    """
    return results_dict