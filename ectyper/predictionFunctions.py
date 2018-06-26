#!/usr/bin/env python

import json
import logging
import os
import pandas as pd

LOG = logging.getLogger(__name__)

"""
    Serotype prediction for E. coli
"""


def predict_serotype(blast_output_file, ectyper_dict_file):
    """
    Predict the serotype of all genomes, given the blast output of the markers against the genomes

    :param blast_output_file: Results of allele file against the genomes of interest
    :param ectyper_dict_file: JSON file of known alleles and their O and H mappings
    :return: The CSV formatted predictions file
    """

    LOG.info("Predicting serotype from blast output")
    output_df = blast_output_to_df(blast_output_file)
    ectyper_df = ectyper_dict_to_df(ectyper_dict_file)
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        LOG.debug("blast_df:\n{}".format(output_df))
        LOG.debug("ectyper_df:\n{}".format(ectyper_df))

    # Merge output_df and ectyper_df
    output_df = output_df.merge(ectyper_df, left_on='qseqid', right_on='name', how='left')
    predictions_dict = {}

    # Select individual genomes
    output_df['genome_name'] = output_df['sseqid'].str.split('|').str[1]

    # Make prediction for each genome based on blast output
    for genome_name, per_genome_df in output_df.groupby('genome_name'):
        predictions_dict[genome_name] = get_prediction(per_genome_df)


    LOG.info("Serotype prediction completed")
    LOG.debug("Predictions dict:\n{}".format(predictions_dict))
    return predictions_dict


def get_prediction(per_genome_df):
    """
     Make serotype prediction for a single genome based on the blast output

    :param per_genome_df: The blastn results for the given genome
    :return: serotype dictionary
    """

    per_genome_df = per_genome_df.sort_values('score', ascending=False)
    LOG.debug("per_genome_df:\n{}".format(per_genome_df))

    # The DataFrame is sorted in descending order by score
    # Once we hit a wzx / wzy or wzm / wzt pair, O type is done
    # Once we hit a flX gene, H type is done

    serotype = {
        'O':'-',
        'H':'-'
    }

    for row in per_genome_df.itertuples():
        # if O or H is already set, skip
        # get the 'O' or 'H' from the antigen column
        ant = row.antigen[:1]
        if serotype[ant] == '-':
            if ant == 'H':
                serotype[ant] = row.antigen
                serotype[row.antigen] = {
                    row.gene:row.score
                }
            else:
                # logic for O-type pairs
                # skip if an allele for a gene already exists
                if row.antigen in serotype and row.gene in serotype[row.antigen]:
                    continue
                else:
                    # if antigen has never been encountered, init
                    if row.antigen not in serotype:
                        serotype[row.antigen] = {}

                    serotype[row.antigen][row.gene] = row.score
                    # if wzm / wzy or wzx / wzy, call the match
                    if 'wzx' in serotype[row.antigen] and 'wzy' in serotype[row.antigen]:
                        serotype[ant] = row.antigen
                    elif 'wzm' in serotype[row.antigen] and 'wzt' in serotype[row.antigen]:
                        serotype[ant] = row.antigen
                    else:
                        continue
        else:
            continue

    return serotype


def blast_output_to_df(blast_output_file):
    """
    Convert the raw Blast output to a DataFrame

    :param blast_output_file: Blast results to convert
    :return: DataFrame of the blast results
    """

    output_data = []
    with open(blast_output_file, 'r') as fh:
        for line in fh:
            fields = line.strip().split()
            entry = {
                'qseqid': fields[0],
                'qlen': fields[1],
                'sseqid': fields[2],
                'length': fields[3],
                'pident': fields[4],
                'sstart': fields[5],
                'send': fields[6],
                'sframe': fields[7],
                'qcovhsp': fields[8],
                'sseq': fields[9]
            }
            output_data.append(entry)
    df = pd.DataFrame(output_data)

    if not output_data:
        LOG.warning("No hits found for blast output file {}".format(blast_output_file))

        # Return empty dataframe with correct columns
        return pd.DataFrame(
            columns=[
                'length', 'pident', 'qcovhsp',
                'qlen', 'qseqid', 'send',
                'sframe', 'sseqid', 'sstart', 'sseq'
            ])
    else:
        df['score'] = df['pident'].astype(float)*df['qcovhsp'].astype(float)/10000
        return df


def ectyper_dict_to_df(ectyper_dict_file):
    """
    Load all the known alleles for the O and H genes, and store them as a DataFrame.
    :param ectyper_dict_file: JSON file of all known O- and H- alleles
    :return: DataFrame of the JSON file
    """

    with open(ectyper_dict_file) as fh:
        ectyper_dict = json.load(fh)
        temp_list = []
        for antigen, alleles in ectyper_dict.items():
            for name, allele in alleles.items():
                new_entry = {
                    'antigen': allele.get('allele'),
                    'name': name,
                    'gene': allele.get('gene'),
                    'desc': allele.get('desc')
                }
                temp_list.append(new_entry)
        df = pd.DataFrame(temp_list)
        return df


def report_result(final_dict, output_file):
    """
    Outputs the results of the ectyper run to the output file, and to the log.

    :param final_dict: Final ectyper predictions dictionary
    :param output_file: File whose contents will be added to the log
    :return: None
    """

    output_line = []
    for k, v in final_dict.items():
        output_line.append(k)

        if 'error' in v:
            output_line.append(v['error'])
        else:
            output_line.append(v['O'])
            output_line.append(v['H'])

            antigens = [v['O'], v['H']]
            for ant in antigens:
                if v[ant] != "-":
                    for kk, vv in sorted(v[ant].items()):
                        output_line.append(kk + ':' + str(vv))

    print_line = "\t".join(output_line)
    with open(output_file, "w") as ofh:
        ofh.write(print_line)
        LOG.info(print_line)


def add_non_predicted(all_genomes_list, predictions_dict, other_dict):
    """
    Add genomes that do not show up in the blast results to the final predictions

    :param all_genomes_list: the list of genomes given by the user
    :param predictions_data_frame: the Dict containing the ectyper predictions
    :return: modified prediction file
    """

    # genome names are given without the filename extension
    for g in all_genomes_list:
        gname = os.path.splitext(os.path.split(g)[1])[0]

        if gname not in predictions_dict:
            LOG.info("gname is {}".format(gname))
            predictions_dict[gname] = {
                'error': other_dict[g]
            }

    return predictions_dict

