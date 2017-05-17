#!/usr/bin/env python

import re
import logging

log = logging.getLogger(__name__)

"""
    Serotype prediction for E. coli
"""


def parse_serotype(blast_record, data):
    """
    Look up the serotype prediction based on the BLAST record.
    The allele should exist in either the O or H database, or there is
    an error, such as "Onovel" being the top match.
    
    
    """

    if blast_record['qseqid'] in data['O']:
        stype = 'O'
    elif  blast_record['qseqid'] in data['H']:
        stype = 'H'
    else:
        return {}

    return {data[stype][blast_record['qseqid']]['gene']: \
                        {'antigen': data[stype][blast_record['qseqid']]['allele'],
                         'blast_record': blast_record}}



def predict_serotype(results_dict):
    """
    Additional logic to get the O and H type from the parsed dictionary of
     blast results. If there are no conflicts between the results, we can use
     them directly.
    :param: results_dict: parsed blast results for serotype prediction 
    :return: results_dict with otype and htype predicted
    """
    log.info("Predicting serotype from parsed results")

    antigen_dict = {'O':re.compile('^(O\d+)'),
                    'H':re.compile('^(H\d+)')}

    for genome_name in results_dict.keys():
        log.debug("Prediction serotype for " + genome_name)

        current_sero_dict = {'O':{'ant_number':None,
                                  'strength':None,
                                  'conflict':{'ant_number':None,
                                              'strength':None},
                                  'gnd':None},
                             'H':{'ant_number':None,
                                  'strength':None,
                                  'conflict':{'ant_number':None,
                                              'strength':None}}}
        if 'serotype' in results_dict[genome_name]:
            for gene_name in results_dict[genome_name]['serotype'].keys():
                current_pident = results_dict[genome_name]['serotype'][gene_name]['blast_record']['pident']

                for antigen in antigen_dict.keys():
                    # get antigen from results
                    m = antigen_dict[antigen].search(results_dict[genome_name]['serotype'][gene_name]['antigen'])
                    if m:
                        match_sero = m.group(1)

                        # We only want to consider 'gnd' if there is conflict or
                        # no other O antigen information
                        if gene_name == 'gnd':
                            current_sero_dict[antigen]['gnd']=match_sero
                            continue

                        if not current_sero_dict[antigen]['ant_number']:
                            current_sero_dict[antigen]['ant_number'] = match_sero
                            current_sero_dict[antigen]['strength']= current_pident
                        elif match_sero == current_sero_dict[antigen]['ant_number']:
                            # good, they agree
                            if current_pident > current_sero_dict[antigen]['ant_number']:
                                current_sero_dict[antigen]['ant_number'] = current_pident
                        else:
                            #bad, they do not agree
                            current_sero_dict[antigen]['conflict']['ant_number'] = match_sero
                            current_sero_dict[antigen]['conflict']['strength'] = current_pident
                            log.debug("Additional prediction required:")
                            log.debug(current_sero_dict)
                    break

            for antigen in current_sero_dict.keys():
                log.debug(current_sero_dict)
                log.debug(antigen)

                if current_sero_dict[antigen]['conflict']['ant_number']:
                    log.debug("Additional serotype prediction required")
                    current_sero_dict = resolve_antigenic_conflict(current_sero_dict)
                    log.debug(current_sero_dict)

                antigen_type = antigen.lower() + "type"

                if current_sero_dict[antigen]['ant_number']:
                    results_dict[genome_name][antigen_type]={'ant_number':current_sero_dict[antigen]['ant_number'],
                                                             'strength':current_sero_dict[antigen]['strength']}
                else:
                    results_dict[genome_name][antigen_type]={'ant_number':'-',
                                                             'strength':0}

    return results_dict


def resolve_antigenic_conflict(sero_dict):
    """
    IF there is conflict, resolve it using additional information if possible,
    :param sero_dict: 
    :return: modified sero_dict, with conflict resolved
    """

    log.debug("Resolving antigenic conflict")


    return sero_dict