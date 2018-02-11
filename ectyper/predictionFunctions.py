#!/usr/bin/env python

from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import absolute_import
from builtins import open
from future import standard_library
standard_library.install_aliases()
import json
import logging
import os
from collections import defaultdict

import pandas as pd

LOG = logging.getLogger(__name__)

"""
    Serotype prediction for E. coli
"""

def df_object_to_unicode(object_df):
    types = object_df.apply(lambda x: pd.api.types.infer_dtype(x.values))
    for col in types[types=='object'].index:
        LOG.warn("Converting {0} object to unicode.".format(col))
        object_df[col] = object_df[col].astype(unicode)
    print(object_df.dtypes)

def predict_serotype(blast_output_file, ectyper_dict_file, predictions_file, detailed=False):
    """Make serotype prediction for all genomes based on blast output

    Args:
        blast_output_file(str):
            blastn output with outfmt:
                "6 qseqid qlen sseqid length pident sstart send sframe qcovhsp -word_size 11"
    ectyper_dict_file(str):
        mapping file used to associate allele id to allele informations
    predictions_file(str):
        csv file to store result
    detailed(bool, optional):
        whether to generate detailed output or not

    Returns:
        predictions_file
    """
    basename, extension = os.path.splitext(predictions_file)
    parsed_output_file = ''.join([basename, '_raw', extension])

    LOG.info("Predicting serotype from blast output")
    output_df = blast_output_to_df(blast_output_file)
    ectyper_df = ectyper_dict_to_df(ectyper_dict_file)
    # Merge output_df and ectyper_df
    output_df = output_df.merge(ectyper_df, left_on='qseqid', right_on='name', how='left')
    predictions_dict = {}
    # Select individual genomes
    output_df['genome_name'] = output_df['sseqid'].str.split('|').str[1]
    # Initialize constants
    gene_pairs = {'wzx':'wzy', 'wzy':'wzx', 'wzm':'wzt', 'wzt':'wzm'}
    predictions_columns = ['O_prediction', 'O_info', 'H_prediction', 'H_info']
    gene_list = ['wzx', 'wzy', 'wzm', 'wzt', 'fliC', 'fllA', 'flkA', 'flmA', 'flnA']
    if detailed:
        # Add gene lists for detailed output report
        for gene in gene_list:
            predictions_columns.append(gene)
    for genome_name, per_genome_df in output_df.groupby('genome_name'):
        # Make prediction for each genome based on blast output
        predictions_dict[genome_name] = get_prediction(
            per_genome_df, predictions_columns, gene_pairs, detailed)
    predictions_df = pd.DataFrame(predictions_dict).transpose()
    if predictions_df.empty:
        predictions_df = pd.DataFrame(columns=predictions_columns)
    predictions_df = predictions_df[predictions_columns]
    # TODO: rm This
    print(output_df.dtypes)
    df_object_to_unicode(output_df)
    print(output_df.dtypes)
    store_df(output_df, parsed_output_file)
    print(predictions_df.dtypes)
    store_df(predictions_df, predictions_file)
    LOG.info("Serotype prediction completed")
    return predictions_file

def get_prediction(per_genome_df, predictions_columns, gene_pairs, detailed, ):
    """Make serotype prediction for single genomes based on blast output

    Args:
        per_genome_df(DataFrame):
            blastn outputs for the given genome
        predictions_columns(dict):
            columns to be filled by prediction function
        gene_pairs(dict):
            dict of pair genes used for paired logic
        detailed(bool):
            whether to generate detailed output or not
    Return:
        Prediction dictionary for the given genome
    """
    # Extract potential predictors
    useful_columns = [
        'gene', 'serotype', 'score', 'name', 'desc', 'pident', 'qcovhsp', 'qseqid', 'sseqid'
    ]
    per_genome_df = per_genome_df.sort_values(['gene', 'serotype', 'score'], ascending=False)
    per_genome_df = per_genome_df[~per_genome_df.duplicated(['gene', 'serotype'])]
    predictors_df = per_genome_df[useful_columns]
    predictors_df = predictors_df.sort_values('score', ascending=False)
    predictions = {}
    for column in predictions_columns:
        predictions[column] = None
    for predicting_antigen in ['O', 'H']:
        genes_pool = defaultdict(list)
        for index, row in predictors_df.iterrows():
            gene = row['gene']
            if detailed:
                predictions[gene] = True
            if not predictions[predicting_antigen+'_prediction']:
                serotype = row['serotype']
                if serotype[0] is not predicting_antigen:
                    continue
                genes_pool[gene].append(serotype)
                prediction = None
                if len(serotype) < 1:
                    continue
                antigen = serotype[0].upper()
                if antigen != predicting_antigen:
                    continue
                if gene in gene_pairs.keys():
                    predictions[antigen+'_info'] = 'Only unpaired alignments found'
                    # Pair gene logic
                    potential_pairs = genes_pool.get(gene_pairs.get(gene))
                    if potential_pairs is None:
                        continue
                    if serotype in potential_pairs:
                        prediction = serotype
                else:
                    # Normal logic
                    prediction = serotype
                if prediction is None:
                    continue
                predictions[antigen+'_info'] = 'Alignment found'
                predictions[predicting_antigen+'_prediction'] = prediction
        # No alignment found
        ## Make prediction based on non-paired gene
        ##   if only one non-paired gene is avaliable
        if predictions.get(predicting_antigen+'_prediction') is None:
            if len(genes_pool) == 1:
                serotypes = list(genes_pool.values())[0]
                if len(serotypes) == 1:
                    predictions[antigen+'_info'] = 'Lone unpaired alignment found'
                    predictions[predicting_antigen+'_prediction'] = serotypes[0]
    return predictions

def blast_output_to_df(blast_output_file):
    '''Convert raw blast output file to DataFrame
    Args:
        blast_output_file(str): location of blast output

    Returns:
        DataFrame:
            DataFrame that contains all informations from blast output
    '''
    # Load blast output
    output_data = []
    with open(blast_output_file, 'r') as fh:
        for line in fh:
            if type(line) is not unicode:
                LOG.warn("Non-unicode line detected in: {0}.".format(blast_output_file))
                line = unicode(line)
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
    if not output_data:
        LOG.info("No hit found for this blast query")
        # Return empty dataframe with correct columns
        return pd.DataFrame(
            columns=[
                'length', 'pident', 'qcovhsp',
                'qlen', 'qseqid', 'send',
                'sframe', 'sseqid', 'sstart'
            ])
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
                for key in new_entry:
                    if type(key) is not unicode:
                        LOG.warn("{0} is of type {1} not unicode.".format(key, type(key)))
                        new_entry[key] = unicode(new_entry[key])
                temp_list.append(new_entry)
        df = pd.DataFrame(temp_list)
        return df

def store_df(src_df, dst_file):
    """Append dataframe to a file if it exists, otherwise, make a new file
    Args:
        src_df(str): dataframe object to be stored
        dst_file(str): dst_file to be modified/created
    """
    if os.path.isfile(dst_file):
        with open(dst_file, 'a') as fh:
            src_df.to_csv(fh, header=False, encoding='utf-8')
    else:
        with open(dst_file, 'w') as fh:
            src_df.to_csv(fh, header=True, index_label='genome', encoding='utf-8')

def report_result(csv_file):
    '''Report the content of dataframe in log

    Args:
        csv_file(str): location of the prediction file
    '''
    df = pd.read_csv(csv_file, encoding='utf-8')
    if df.empty:
        LOG.info('No prediction was made becuase no alignment was found')
        return
    LOG.info('\n{0}'.format(df.to_string(index=False)))

def add_non_predicted(all_genomes_list, predictions_file):
    '''Add genomes that do not show up in blast result to prediction file

    Args:
        all_genome_list(list):
            list of genomes from user input
        predictions_file(str):
            location of the prediction file
    Returns:
        str: location of the prediction file
    '''
    df = pd.read_csv(predictions_file, encoding='utf-8')
    df = df.merge(pd.DataFrame(all_genomes_list, columns=['genome']), on='genome', how='outer')
    df.fillna({'O_info':'No alignment found', 'H_info':'No alignment found'}, inplace=True)
    df.fillna('-', inplace=True)
    df.to_csv(predictions_file, index=False, encoding='utf-8')
    return predictions_file
