#!/usr/bin/env python

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
    elif blast_record['qseqid'] in data['H']:
        stype = 'H'
    else:
        log.warning("{0} has no match in either O or H database".format(
            blast_record['qseqid']))
        return {}

    return {data[stype][blast_record['qseqid']]['gene']: \
                {'antigen': data[stype][blast_record['qseqid']]['allele'],
                 'blast_record': blast_record,
                 'stype': stype}}


def predict_serotype(results_dict):
    """
    Additional logic to get the O and H type from the parsed dictionary of
     blast results. If there are no conflicts between the results, we can use
     them directly.
    :param: results_dict: parsed blast results for serotype prediction
    :return: results_dict with otype and htype predicted
    """
    log.info("Predicting serotype from parsed results")

    for genome_name in results_dict.keys():
        log.debug("Predicting serotype for " + genome_name)

        sero_dict = {'O': '-', 'H': '-'}
        for gene_name in results_dict[genome_name]['serotype']:
            # skip if gnd, we only want it in cases of disagreement
            if gene_name == 'gnd':
                pass
            else:
                stype = results_dict[genome_name]['serotype'][gene_name][
                    'stype']
                antigen = results_dict[genome_name]['serotype'][gene_name][
                    'antigen']

                if stype == 'O':
                    pass

                if sero_dict[stype] == '-':
                    sero_dict[stype] = antigen
                else:
                    if sero_dict[stype] == antigen:
                        # good, as expected
                        pass
                    else:
                        if gene_name == 'gnd':
                            # last resort to use gnd, only if other conflicts
                            pass
                        else:
                            log.warning("{0}type for {1} conflicts! {2} {3}"
                                        "".format(stype, genome_name,
                                                  sero_dict[stype], antigen))

        results_dict[genome_name]['serotype']['otype'] = sero_dict['O']
        results_dict[genome_name]['serotype']['htype'] = sero_dict['H']
    
        if sero_dict['O'] == sero_dict['H'] == '-':
            log.warning("No serotype found for %s!!!", genome_name)
    return results_dict


def resolve_antigenic_conflict(sero_dict):
    """
    IF there is conflict, resolve it using additional information if possible,
    :param sero_dict:
    :return: modified sero_dict, with conflict resolved
    """

    log.debug("Resolving antigenic conflict")

    return sero_dict
