#!/usr/bin/env python

import json
import logging
import os
import pandas as pd

LOG = logging.getLogger(__name__)

"""
    Serotype prediction for E. coli
"""


def predict_serotype(blast_output_file, ectyper_dict_file, args):
    """
    Predict the serotype of all genomes, given the blast output of the markers against the genomes

    :param blast_output_file: Results of allele file against the genomes of interest
    :param ectyper_dict_file: JSON file of known alleles and their O and H mappings
    :param args: Commandline arguments
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
        predictions_dict[genome_name] = get_prediction(per_genome_df, args)


    LOG.info("Serotype prediction completed")
    LOG.debug("Predictions dict:\n{}".format(predictions_dict))
    return predictions_dict


def get_prediction(per_genome_df, args):
    """
     Make serotype prediction for a single genome based on the blast output

    :param per_genome_df: The blastn results for the given genome
    :param args: Commandline args
    :return: serotype dictionary
    """

    per_genome_df = per_genome_df.sort_values(by=['score'], ascending=False)
    LOG.debug("per_genome_df:\n{}".format(per_genome_df))

    # The DataFrame is sorted in descending order by score
    # Once we hit a wzx / wzy or wzm / wzt pair, O type is done
    # Once we hit a flX gene, H type is done

    serotype = {
        'O':'-',
        'H':'-'
    }

    otype = {}
    # Go for the highest match, if both genes exist over the thresholds
    best_order = []
    for row in per_genome_df.itertuples():
        # H is already set, skip
        # get the 'O' or 'H' from the antigen column
        ant = row.antigen[:1]

        if ant == 'H' and serotype[ant] == '-':
            serotype[ant] = row.antigen
            serotype[row.antigen] = {
                row.gene:row.score
            }
            if args.sequence:
                serotype[row.antigen]["≈" + row.gene] = row.sseq
        elif ant == 'O':
            # logic for O-type pairs
            # skip if an allele for a gene already exists
            if row.antigen in otype and row.gene in otype[row.antigen]:
                continue
            else:
                # if antigen has never been encountered, init
                if row.antigen not in otype:
                    otype[row.antigen] = {}
                    best_order.append(row.antigen)

                if args.sequence:
                    otype[row.antigen]["≈" + row.gene] = row.sseq

                otype[row.antigen][row.gene] = row.score


    # having gone through all the hits over the threshold, make the call
    # go through the O-antigens in order, making the call on the first that have
    # a matching pair
    for o in best_order:
        # if wzm / wzy or wzx / wzy, call the match
        if 'wzx' in otype[o] and 'wzy' in otype[o]:
            serotype['O'] =  o
            serotype[o]=otype[o]
            break
        elif 'wzm' in otype[o] and 'wzt' in otype[o]:
            serotype['O'] = o
            serotype[o] = otype[o]
            break

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

    output = []
    for k, v in final_dict.items():
        output_line = [k]

        if 'error' in v:
            output_line.append(v['error'])
        else:
            output_line.append(v['O'])
            output_line.append(v['H'])

            antigens = [v['O'], v['H']]
            for ant in antigens:
                if ant != "-":
                    for kk, vv in sorted(v[ant].items()):
                        if "≈" in kk:
                            output_line.append(kk + ':' + vv)
                        else:
                            output_line.append(kk + ':' + " {0:.2f}".format(vv))

        print_line = "\t".join(output_line)
        output.append(print_line + "\n")
        LOG.info(print_line)

    with open(output_file, "w") as ofh:
        for line in sorted(output):
            ofh.write(line)


def add_non_predicted(all_genomes_list, predictions_dict, other_dict):
    """
    Add genomes that do not show up in the blast results to the final predictions

    :param all_genomes_list: the list of genomes given by the user
    :param predictions_data_frame: the Dict containing the ectyper predictions
    :return: modified prediction file
    """

    # test on '/mnt/moria/enterobase_serotype/ESC_GA9165AA_AS.fasta'
    # genome names are given without the filename extension
    for g in all_genomes_list:
        gname = os.path.splitext(os.path.split(g)[1])[0]

        if gname not in predictions_dict:
            if g in other_dict:
                predictions_dict[gname] = {
                    'error': other_dict[g]
                }
            else:
                predictions_dict[gname] = {
                    'error': "No serotyping-specific genes found"
                }

    return predictions_dict
