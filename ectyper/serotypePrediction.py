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

        sero_dict = {'O': '-', 'H': '-', 'O allele': '-', 'H allele': '-'}

        # vars for paired gene logic
        gene_pairs = {'wzx':'wzy', 'wzy':'wzx', 'wzm':'wzt', 'wzt':'wzm'}
        existing_gene_dict = {}

        for gene_name in results_dict[genome_name]['serotype']:

            # skip if gnd, we only want it in cases of disagreement
            if gene_name == 'gnd':
                continue

            # Obtain information of the best allele for this gene
            stype = results_dict[genome_name]['serotype'][gene_name][
                'stype']
            antigen = results_dict[genome_name]['serotype'][gene_name][
                'antigen']
            blast_record = results_dict[genome_name]['serotype'][gene_name][
                'blast_record']
            allele = blast_record['qseqid']

            score = (float(blast_record['pident'])/100) * (float(blast_record['length'])/float(blast_record['qlen']))

            # Special logic for paired gene
            # get the name of the other gene
            other_gene_name = gene_pairs.get(gene_name)
            if other_gene_name:
                # add this gene to existing_gene_dict if it is one of the gene pairs
                if other_gene_name:
                    entry = {
                        'stype': stype,
                        'antigen': antigen,
                        'allele': allele,
                        'score': score
                    }
                    existing_gene_dict[gene_name]=entry
                # check if the other gene is already in the dictionary
                other_gene_allele = existing_gene_dict.get(other_gene_name)
                if other_gene_allele is None:
                    continue
                # add the allele only if the stype and antigen matches
                if (other_gene_allele.get('stype') == stype) & \
                (other_gene_allele.get('antigen') == antigen):
                    pass
                else:
                    continue
            
            # use this allele if there is no allele yet
            if sero_dict[stype] == '-':
                sero_dict[stype] = antigen
                sero_dict[stype+' allele'] = allele
                continue
            if sero_dict[stype] == antigen:
                # good, as expected
                continue
            # conflict between alleles
            log.warning("{0}type for {1} conflicts! {2} {3}"
                        "".format(stype, genome_name,
                                    sero_dict[stype], antigen))

        # Logic when only one of the paired allele is found
        for stype in ['O', 'H']:
            cur_type = sero_dict[stype]
            if cur_type == '-':
                # Use the best non-paired pair gene if nothing else can be used
                unpaired_alleles = existing_gene_dict.values()
                if len(unpaired_alleles) > 0:
                    # Choose the one with highest (percent identitity * query coverage)
                    sorted_alleles = sorted(unpaired_alleles, key=lambda k: k['score'], reverse=True)
                    target_allele = sorted_alleles[0]
                    sero_dict[stype] = target_allele['antigen']
                    sero_dict[stype+' allele'] = target_allele['allele']

            results_dict[genome_name]['serotype'][stype.lower()+'type'] = sero_dict[stype]
            results_dict[genome_name]['serotype'][stype+' allele'] = sero_dict[stype+' allele']

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
