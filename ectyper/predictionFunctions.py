#!/usr/bin/env python

import logging
import os
import pandas as pd
import json
LOG = logging.getLogger(__name__)

"""
    Serotype prediction for E. coli
"""

def predict_serotype(blast_output_file, ectyper_dict_file, predictions_file):
    """
    Make serotype prediction based on blast output
    :param blast_output_file: blastn output with 
        outfmt "6 qseqid qlen sseqid length pident sstart send sframe qcovhsp -word_size 11"
    :param ectyper_dict_file: mapping file used to associate allele id to allele informations
    :param predictions_file: csv file to store result
    :return: predictions_file
    """
    basename, extension = os.path.splitext(predictions_file)
    parsed_output_file = ''.join([basename, '_raw', extension])
    useful_columns = [
        'gene', 'serotype', 'score', 'name', 'desc', 'pident', 'qcovhsp', 'qseqid', 'sseqid'
    ]

    LOG.info("Predicting serotype from blast output")
    output_df = blast_output_to_df(blast_output_file)
    store_df(output_df, parsed_output_file)
    ectyper_df = ectyper_dict_to_df(ectyper_dict_file)
    # Merge output_df and ectyper_df
    output_df = output_df.merge(ectyper_df, left_on='qseqid', right_on='name', how='left')
    predictions_dict = {}
    # Select individual genomes
    for name, per_genome_df in output_df.groupby('sseqid'):
        df = per_genome_df
        # Extract potential predictors
        df = df.sort_values(['gene', 'serotype', 'score'])
        df = df[~df.duplicated(['gene', 'serotype'])]
        predictors_df = df[useful_columns]
        predictors_df = predictors_df.sort_values('score', ascending=False)
        predictors_df.to_csv('output.csv')

        # Make prediction based on predictors
        from collections import defaultdict
        GENE_PAIRS = {'wzx':'wzy', 'wzy':'wzx', 'wzm':'wzt', 'wzt':'wzm'}
        predictions = {
            'H': None,
            'O': None
        }
        for predicting_antigen in ['O', 'H']:
            genes_pool = defaultdict(list)
            for index, row in predictors_df.iterrows():
                gene = row['gene']
                serotype = row['serotype']
                genes_pool[gene].append(serotype)
                prediction = None
                if len(serotype) < 1:
                    continue
                antigen = serotype[0].upper()
                if antigen != predicting_antigen:
                    continue
                if gene in GENE_PAIRS.keys():
                    # Pair gene logic
                    potential_pairs = genes_pool.get(GENE_PAIRS.get(gene))
                    if potential_pairs is None:
                        continue
                    if serotype in potential_pairs:
                        prediction = serotype
                else:
                    # Normal logic
                    prediction = serotype
                if prediction==None:
                    continue
                predictions[predicting_antigen]=prediction
                break
        predictions_dict[name] = predictions
    predictions_df = pd.DataFrame(predictions_dict).transpose()
    store_df(predictions_df, predictions_file)
    LOG.info("Serotype prediction completed")
    return predictions_file

def blast_output_to_df(blast_output_file):
    # Load blast output
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
                'qcovhsp': fields[8]
            }
            output_data.append(entry)
    df = pd.DataFrame(output_data)
    df['score'] = df['pident'].astype(float)*df['qcovhsp'].astype(float)/10000
    return df

def ectyper_dict_to_df(ectyper_dict_file):
    # Read ectyper dict
    with open(ectyper_dict_file) as fh:
        ectyper_dict = json.load(fh)
        temp_list = []
        for antigen, alleles in ectyper_dict.items():
            for name, allele in alleles.items():
                new_entry = {
                    'serotype': allele.get('allele'),
                    'name': name,
                    'gene': allele.get('gene'),
                    'desc': allele.get('desc')
                }
                temp_list.append(new_entry)
        df = pd.DataFrame(temp_list)
        return df

def store_df(src_df, dst_file):
    """
    Append dataframe to a file if it exists, otherwise, make a new file
    :param src_df: dataframe object to be stored
    :param dst_file: dst_file to be modified/created
    """
    if os.path.isfile(dst_file):
        with open(dst_file, 'a') as fh:
            src_df.to_csv(fh, header=False)
    else:
        with open(dst_file, 'w') as fh:
            src_df.to_csv(fh, header=True, index_label='index')

def report_result(csv_file):
    # Report the content of dataframe in log
    LOG.info('\n%s',pd.read_csv(csv_file).to_string())