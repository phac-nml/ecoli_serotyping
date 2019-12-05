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
    output_df = blast_output_to_df(blast_output_file) #columns: length	pident	qcovhsp	qlen qseqid	send	sframe	sseq	sseqid	sstart	score
    if args.debug:
        LOG.debug("Wrote BLAST results file blast_outputBLAST_df.tsv in {}".format(os.getcwd()))
        output_df.to_csv("blast_outputBLAST_df.tsv",sep="\t") #DEBUG

    ectyper_df = ectyper_dict_to_df(ectyper_dict_file) #columns: antigen	desc	gene	name
    #ectyper_df.to_csv("ectyper_df.tsv", sep="\t")  # DEBUG
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        LOG.debug("blast_df:\n{}".format(output_df))
        LOG.debug("ectyper_df:\n{}".format(ectyper_df))

    # Merge output_df and ectyper_df
    output_df = output_df.merge(ectyper_df, left_on='qseqid', right_on='name', how='left')
    predictions_dict = {}


    # Select individual genomes
    output_df['genome_name'] = output_df['sseqid'].str.split('|').str[1] #e.g. lcl|Escherichia_O26H11|17

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

    per_genome_df = per_genome_df.sort_values(by=['score'], ascending=False,
                                              kind='mergesort')

    if args.debug:
        per_genome_df.to_csv("per_genome_df.tsv", sep="\t")
        LOG.debug("Wrote per_genome_df.tsv to {}".format(os.getcwd()))

    # The DataFrame is sorted in descending order by score e.g. serotype={'O': 'O26', 'H': 'H11', 'H11': {'fliC': 1.0}, 'O26': {'wzx': 0.46}}
    serotype = {
        'O':'-',
        'H':'-'
    }

    otype = {} #{antigen:{gene:score}} key value pairs e.g. otype ={'O26':{'wzx':1.00,'wzy':1.00}} or {'O26': {'wzx': 0.46}}
    # Go for the highest match, if both genes exist over the thresholds

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
                    #best_order.append(row.antigen)
                if args.sequence:
                    otype[row.antigen]["≈" + row.gene] = row.sseq

                otype[row.antigen][row.gene] = row.score

    rank_Otype_dict={}
    #rank O-type serovars based on the sum of scores from BOTH alleles (wzx/wzy or wzt/wzm)
    for otypename in otype.keys():
        rank_Otype_dict[otypename]={"numalleles":0, "scores":[], "sumscore":0}
        alleles = [key for key in otype[otypename].keys() if "≈" not in key]
        rank_Otype_dict[otypename]["numalleles"] = len(alleles)
        scores = [otype[otypename][allele]  for allele in alleles]
        rank_Otype_dict[otypename]["scores"] = scores
        rank_Otype_dict[otypename]["sumscore"] = sum(scores)


    #maxscore=0; best_order_list = [] #based on score centered on identity and coverage e.g. best_order = ['O26', 'O174']
    scorestupleslist = [(otypename,rank_Otype_dict[otypename]["sumscore"]) for otypename in rank_Otype_dict.keys()]
    scorestupleslist = sorted(scorestupleslist, key=lambda x: x[1], reverse=True)
    best_order_list = [item[0] for item in scorestupleslist]


    LOG.debug("Otype dict:\n{}".format(otype))
    LOG.debug("Serotype dict:\n{}".format(serotype))
    LOG.debug("\"Best order alleles-scores\" list:\n{}".format(scorestupleslist))
    LOG.debug("\"Best order alleles\" list:\n{}".format(best_order_list))

    # having gone through all the hits over the threshold, make the call
    # go through the O-antigens in order, making the call on the first that have
    # a matching pair
    for o in best_order_list:
        # if wzm / wzy or wzx / wzy, call the match
        if 'wzx' in otype[o] and 'wzy' in otype[o]:
            serotype['O'] =  o
            serotype[o]=otype[o]
            break
        elif 'wzm' in otype[o] and 'wzt' in otype[o]:
            serotype['O'] = o
            serotype[o] = otype[o]
            break
        # FIX: O-antigen typing might fail due to poor sequencing or inability to assemble one of the wzx/wzy/wzm/wzt loci
        # if only one of the signatures is found, still produce output but warn user on false positives
        elif 'wzx' in otype[o] or 'wzy' in otype[o]:
            serotype['O'] = o
            serotype[o] = otype[o]
            break
        elif 'wzm' in otype[o] or 'wzt' in otype[o]:
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

def mean(numbers):
    return sum(numbers)/len(numbers)

def quality_control_results(sample, final_results_dict):
    """
    Determined approximate quality of the prediction based on the allele scores. Adopt pessimistic approach by looking at min values
    :param sample: sample/genome name)
    :param final_results_dict: dictionary with final output results (e.g. serovar, sequences, conf. scores)
    :return: Quality dictionary
    """

    Otype=final_results_dict[sample]["O"]; Oscores=[]
    Htype=final_results_dict[sample]["H"]; Hscore  = 0
    AllelsList=[]
    species=final_results_dict[sample]["species"]

    if Otype != "-":
        for allele in final_results_dict[sample][Otype].keys():
            if any([item == allele for item in ["wzx","wzy","wzm","wzt"]]):
              Oscores.append(final_results_dict[sample][Otype][allele])
              AllelsList.append(allele)

    if Htype != "-":
        for allele in final_results_dict[sample][Htype].keys():
            if any([item == allele for item in ["fliC"]]):
              Hscore = final_results_dict[sample][Htype][allele]
              AllelsList.append(allele)

    if species != "Escherichia coli":
        return {"QCflag":"NA","AlleleNames":AllelsList,"NumberOfAlleles":len(AllelsList),"ConfidenceLevel":"-"}

    if Otype != "-" and Htype != "-":
        scores = [mean(Oscores),Hscore]
        aggregatescore = mean(scores)
    elif Otype != "-" and Htype == "-":
        aggregatescore = min(Oscores)
    else:
        return {"QCflag":"FAIL","AlleleNames":AllelsList,"NumberOfAlleles":len(AllelsList),"ConfidenceLevel":"-"} #if O and H serovar is not determined OR O is not determined = automatic fail


    if aggregatescore >= 0.95:
        qcflag = "PASS"
        confidencelevel="HIGH"
    elif 0.80 <= aggregatescore < 0.95:
        qcflag = "PASS"
        confidencelevel = "MEDIUM"
    elif 0.50 <= aggregatescore < 0.80:
        qcflag = "PASS"
        confidencelevel = "LOW"
    else:
        qcflag = "FAIL"
        confidencelevel = "-"

    return {"QCflag": qcflag, "AlleleNames": AllelsList, "NumberOfAlleles": len(AllelsList), "ConfidenceLevel":confidencelevel}


def report_result(final_dict, output_dir, output_file):
    """
    Outputs the results of the ectyper run to the output file, and to the log.

    :param final_dict: Final ectyper predictions dictionary {'Sample': {'O': 'O26', 'H': 'H11', 'H11': {'fliC': 1.0}, 'O26': {'wzx': 0.46}}}
    :param output_file: File whose contents will be added to the log
    :return: None
    """

    header = "Name\tSpecies\tO-type\tH-type\tSerovar\tECtyperQC\tConfidence\tEvidence\tAlleles\tWarnings\n"
    output = []
    LOG.info(header)

    alleleseqs={}

    #print(final_dict)
    for sample in final_dict.keys():
        output_line = [sample] #name of a query sample/genome
        output_line.append(final_dict[sample]["species"]) #add species info
        Otype="-"; Htype="-"

        if "O" in final_dict[sample].keys():
            Otype=final_dict[sample]["O"]
        output_line.append(Otype)

        if "H" in final_dict[sample].keys():
            Htype=final_dict[sample]["H"]
        output_line.append(Htype)
        output_line.append("{}:{}".format(Otype,Htype))

        allelesscoresstr = ""
        if Otype != "-":
            for Oallele in sorted(final_dict[sample][Otype].keys()): #makes sure that wzx/wzy and wzm/wzt are reported alphabetically
                if "≈" not in Oallele:
                    score=final_dict[sample][Otype][Oallele]
                    allelesscoresstr = allelesscoresstr + "{}:{:.3f};".format(Oallele, score)
                else:
                    alleleseqs[Oallele+"-"+Otype]=final_dict[sample][Otype][Oallele]

        if Htype != "-":
            for Hallele, score in final_dict[sample][Htype].items():
                if "≈" not in Hallele:
                    allelesscoresstr = allelesscoresstr + "{}:{:.3f};".format(Hallele, score)
                else:
                    alleleseqs[Hallele+"-"+Htype]=final_dict[sample][Htype][Hallele]

        if allelesscoresstr == "":
            allelesscoresstr = "-"
            output_line = output_line + ["-"]*3
        else:
            QCdict = quality_control_results(sample, final_dict)
            output_line.append(QCdict["QCflag"]) #QC flag
            output_line.append(QCdict["ConfidenceLevel"])  #Confidence level
            output_line.append("Based on {} allele(s)".format(QCdict["NumberOfAlleles"])) #evidence
        output_line.append(allelesscoresstr) #allele markers with the corresponding confidence score ranging from 0 to 1
        output_line.append(final_dict[sample]["error"])


        if alleleseqs:
            with open(file=output_dir+"/queryalleleseqs.fasta", mode="w") as fp:
                for allele, seq in alleleseqs.items():
                    fp.write(">{}\n{}\n".format(allele[1:],seq))
                fp.close()

        #if 'error' in v:
        #    output_line.append(v['error'])
        #else:
        #    output_line.append(v['O'])
        #    output_line.append(v['H'])
        #    output_line.append("PASS") #QC flag

        #    antigens = [v['O'], v['H']]
        #    for ant in antigens:
        #        if ant != "-":
        #            for kk, vv in sorted(v[ant].items()):
        #                #if "≈" in kk:
                        #    output_line.append(kk + ':' + vv)
                        #else:
        #               output_line.append(kk + ':' + " {0:.2f}".format(vv))
        #       else:
        #           output_line.append(["-"]*4)

        print_line = "\t".join(output_line)
        output.append(print_line + "\n")
        LOG.info(print_line)

    with open(output_file, mode="w", encoding="utf8") as ofh:
        ofh.write(header)
        for line in sorted(output):
            ofh.write(line)


def add_non_predicted(all_genomes_list, predictions_dict, other_dict, ecoli_dict):
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
            if gname in other_dict:
                predictions_dict[gname] = {
                    'error': other_dict[gname]["error"],
                    'species': other_dict[gname]["species"]
                }
            else:
                predictions_dict[gname] = {
                    'error':  "No serotyping-specific E.coli genes found",
                    'species': ecoli_dict[gname]["species"]
                }

    return predictions_dict
