#!/usr/bin/env python

import re
import src.blastFunctions
import logging

log = logging.getLogger(__name__)

"""
    Serotype prediction for E. coli
"""


def parse_serotype(blast_record, args):
    """
    Entry point for serotype prediction

    :param blast_record: = {'qseqid':la[0],
                        'qlen':la[1],
                        'sseqid':la[2],
                        'length':la[3],
                        'pident':la[4]
                        }
    :param args: argparse commandline options
    :return: dict = {
        otype = O
        htype = H
    }
    """


    # serotype header formats
    # >1__fliC__fliC-H1__1 AB028471.1;flagellin;H1
    # >2__flkA__flkA-H3__77 AB128916;flagellin;H3
    # >3__fllA__fllA-H44__82 #VALUE!;flagellin;H44
    # >4__flmA__flmA-H54__84 AB128918;flagellin;H54
    # >5__flnA__flnA-H17__85 CP002291;flagellin;H17
    # >6__wzm__wzm-O101-Gp15__86 AB812046.1;O-antigen ABC transporter permease;O101
    # >7__wzt__wzt-O101-Gp15__99 AB812046.1;ATP-binding protein;O101
    # >8__wzx__wzx-O1__112 GU299791.1;O antigen flippase;O1
    # >9__wzy__wzy-O1__314 GU299791.1;O antigen polyermase;O1
    # >8__wzx__wzx-O116var1__513 ERS085345;O antigen flippase;O116-like
    # >6__wzm__wzm-Onovel32__523 SRR2544773;O-antigen ABC transporter permease;Onovel32
    # >7__wzt__wzt-Onovel32__524 SRR2544773;ATP-binding protein;Onovel32
    # >fliC_54_AJ605766_H17
    # >wzx_36_DQ462205-FJ539194_O28ac
    # >gnd|1_O26
    #

    serotype_results = {}

    # Initially check that the result passes the length / identity filters
    if not src.blastFunctions.record_passes_cutoffs(blast_record, args):
        return serotype_results

    h_re_patterns = (
        # Look for any flagella match
        re.compile('(fl\w\w)-(H\d+)__'),
        re.compile('(fl\w\w)_.+_(H\d+)$')
    )

    o_re_patterns = (
        # look for any somatic match
        re.compile('(wz\w)-(O\d+)'),
        re.compile('^(wz\w)_.+_(O\d+)\w*$'),
        re.compile('(gnd)\|\d+_(O\d+)')
    )

    all_re_patters = (h_re_patterns, o_re_patterns)

    # Once a match has been found, we don't need to check any others
    # Checks for H-antigens first
    for pattern in all_re_patters:
        for rep in pattern:
            m = rep.search(blast_record['qseqid'])

            if m:
                serotype_results['serotype'] = \
                    {m.group(1): {'antigen':m.group(2),
                                  'blast_record':blast_record}}

                log.debug(serotype_results)
                return serotype_results

    # If no matches, return blank dictionary
    return serotype_results



def predict_serotype(results_dict):
    """
    Additional logic to get the O and H type from the parsed dictionary of
     blast results. If there are no conflicts between the results, we can use
     them directly.
    :param: results_dict: parsed blast results for serotype prediction 
    :return:  
    """
    log.info("Predicting serotype from parsed results")

    antigen_dict = {'O':re.compile('^(O\d+)'),
                    'H':re.compile('^(H\d+)')}

    for genome_name in results_dict.keys():
        log.debug("Prediction serotype for " + genome_name)

        current_sero_dict = {'O':None,
                             'H':None,
                             'additional_prediction':None,
                             'match':None}

        for gene_name in results_dict[genome_name]['serotype'].keys():
            # We only want to consider 'gnd' if there is conflict or
            # no other O antigen information
            if gene_name == 'gnd':
                continue

            current_match = results_dict[genome_name]['serotype'][gene_name]['blast_record']['pident']

            for antigen in antigen_dict.keys():
                # get antigen from results
                m = antigen_dict[antigen].search(results_dict[genome_name]['serotype'][gene_name]['antigen'])
                if m:
                    match_sero = m.group(1)
                    if not current_sero_dict[antigen]:
                        current_sero_dict[antigen] = m.group(1)
                        current_sero_dict['match']= current_match
                    elif match_sero == current_sero_dict[antigen]:
                        # good, they agree
                        # update the match if better
                        if current_match > current_sero_dict['match']:
                            current_sero_dict['match'] = current_match
                    else:
                        #bad, they do not agree
                        current_sero_dict['additional_prediction']=True
                        current_sero_dict[antigen] = current_sero_dict[antigen] + "," + match_sero
                        log.debug("Additional prediction required:")
                        log.debug(match_sero)
                        log.debug(current_sero_dict)
                break

        if current_sero_dict['additional_prediction']:
            log.debug("Additional serotype prediction")
            log.debug(current_sero_dict)
        else:
            if current_sero_dict['O']:
                results_dict[genome_name]['otype']=current_sero_dict['O']
            else:
                results_dict[genome_name]['otype']='NA'

            if current_sero_dict['H']:
                results_dict[genome_name]['htype']=current_sero_dict['H']
            else:
                results_dict[genome_name]['htype']='NA'

    return results_dict
