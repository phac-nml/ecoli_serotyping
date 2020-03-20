#!/usr/bin/env python

import json
import logging
import os, re
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
        LOG.debug("Wrote BLAST results output file in {}".format(os.getcwd()))
        output_df.to_csv("blast_outputBLAST_df.tsv",sep="\t") #DEBUG

    ectyper_df = ectyper_dict_to_df(ectyper_dict_file) #columns: antigen	desc	gene	name

    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    #    LOG.debug("blast_df:\n{}".format(output_df))
    #    LOG.debug("ectyper_df:\n{}".format(ectyper_df))

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

    return predictions_dict, output_df


def get_prediction(per_genome_df, args):
    """
     Make serotype prediction for a single genome based on the blast output

    :param per_genome_df: The blastn reference allele hits results for a given genome
    :param args: Commandline args
    :return: serotype dictionary
    """

    # The per_genome_df DataFrame is sorted in allele score descending order
    per_genome_df = per_genome_df.sort_values(by=['score'], ascending=False,
                                              kind='mergesort')

    serotype_dictkeys = {"serogroup","genescores"}

    serotype = {
        'O':dict.fromkeys(serotype_dictkeys,"-"),
        'H':dict.fromkeys(serotype_dictkeys,"-")
    }

    serotype["O"]["genescores"] = {}; serotype["O"]["alleles"] = {}
    serotype["H"]["genescores"] = {}; serotype["H"]["alleles"] = {}
    serotype["gene2allele"] = {}


    # Go for the highest match, if both genes exist over the thresholds
    blastresultsdict={}; blastresultsdict["O"]={}; blastresultsdict["H"]={}

    for row in per_genome_df.itertuples():
        # H is already set, skip
        # get the 'O' or 'H' from the antigen column
        ant = row.antigen[:1]

        blastresultsdict[ant][row.qseqid] = {}
        blastresultsdict[ant][row.qseqid]["gene"] = row.gene
        blastresultsdict[ant][row.qseqid]["antigen"] = row.antigen
        blastresultsdict[ant][row.qseqid]["score"] = row.score
        blastresultsdict[ant][row.qseqid]["identity"] = row.pident
        blastresultsdict[ant][row.qseqid]["coverage"] = row.qcovhsp
        blastresultsdict[ant][row.qseqid]["contigname"] = row.sseqid.split('|')[2]
        blastresultsdict[ant][row.qseqid]["startpos"] = row.sstart
        blastresultsdict[ant][row.qseqid]["endpos"] = row.send
        blastresultsdict[ant][row.qseqid]["length"] = row.length


    allelefieldnames = ["identity", "coverage", "contigname", "length", "startpos", "endpos", "gene"]

    ant = "H"
    if blastresultsdict[ant].keys():
        topHallele = max([(dballele,blastresultsdict[ant][dballele]["score"]) for dballele in blastresultsdict["H"].keys()])[0]
        serotype[ant]["serogroup"] = blastresultsdict[ant][topHallele]["antigen"]
        gene = blastresultsdict[ant][topHallele]["gene"]
        score = blastresultsdict[ant][topHallele]["score"]
        serotype[ant]["genescores"] = {gene:score}
        serotype["gene2allele"][gene] = topHallele
        serotype[ant]["alleles"][topHallele] = setAlleleMeta(allelefieldnames,blastresultsdict[ant][topHallele])
    else:
        LOG.warning("No H antigen alleles were identified in this sample")

    sortedOalleles = [tuple[0] for tuple in sorted([(dballele,
                                         blastresultsdict["O"][dballele]["score"])
                                         for dballele in blastresultsdict["O"].keys()],
                                         key=lambda x: x[1], reverse=True)]

    otype={}
    for allele in sortedOalleles:
        oantigen = blastresultsdict["O"][allele]["antigen"]
        if blastresultsdict["O"][allele]["antigen"] not in otype.keys():
            otype[oantigen] = {"genescores":{}, "alleles":[], "allele2gene":{}}
        if blastresultsdict["O"][allele]["gene"] not in otype[oantigen]["genescores"].keys():
            gene = blastresultsdict["O"][allele]["gene"]
            otype[oantigen]["genescores"][gene] = blastresultsdict["O"][allele]["score"]
            otype[oantigen]["alleles"].append(allele)
            otype[oantigen]["allele2gene"][allele] = gene


    #rank O-type serovars based on the sum of scores from BOTH alleles (wzx/wzy or wzt/wzm)
    rank_Otype_dict = {}
    for oantigen in otype.keys():
        rank_Otype_dict[oantigen]={"numalleles":0, "scores":[], "sumscore":0}
        rank_Otype_dict[oantigen]["numalleles"] = len(otype[oantigen]["alleles"])
        rank_Otype_dict[oantigen]["scores"] = otype[oantigen]["genescores"].values()
        rank_Otype_dict[oantigen]["sumscore"] = sum(rank_Otype_dict[oantigen]["scores"] )


    scorestupleslist = [(otypename,rank_Otype_dict[otypename]["sumscore"]) for otypename in rank_Otype_dict.keys()]
    scorestupleslist = sorted(scorestupleslist, key=lambda x: x[1], reverse=True) #[('O102', 1.73)]
    best_order_list = [item[0] for item in scorestupleslist] #['O102']


    LOG.debug("Otype dict:\n{}".format(otype))
    LOG.debug("Serotype dict:\n{}".format(serotype))
    LOG.debug("\"Best order alleles-scores\" list:\n{}".format(scorestupleslist))
    LOG.debug("\"Best order alleles\" list:\n{}".format(best_order_list))

    # having gone through all the hits over the threshold, make the call
    # go through the O-antigens in order, making the call on the first that have
    # a matching pair

    ant="O"; selectedOantigen = ""
    for oantigen in best_order_list:
        # if wzm / wzy or wzx / wzy, call the match
        if 'wzx' in otype[oantigen]["genescores"].keys() and 'wzy' in otype[oantigen]["genescores"].keys():
            serotype[ant]['serogroup'] = oantigen
            serotype[ant]['genescores'] = {'wzx':otype[oantigen]["genescores"]["wzx"],
                                           'wzy':otype[oantigen]["genescores"]["wzy"]}
            selectedOantigen = oantigen
            break
        elif 'wzm' in otype[oantigen]["genescores"].keys() and 'wzt' in otype[oantigen]["genescores"].keys():
            serotype[ant]['serogroup'] = oantigen
            serotype[ant]['genescores'] = {'wzm': otype[oantigen]["genescores"]["wzm"],
                                           'wzt': otype[oantigen]["genescores"]["wzt"]}
            selectedOantigen = oantigen
            break

        # FIX: O-antigen typing might fail due to poor sequencing or inability to assemble one of the wzx/wzy/wzm/wzt loci
        # if only one of the signatures is found, still produce output but warn user on false positives
        elif 'wzx' in otype[oantigen]["genescores"].keys() or 'wzy' in otype[oantigen]["genescores"].keys():
            serotype['O']['serogroup']  = oantigen

            if 'wzx' in otype[oantigen]["genescores"].keys():
                serotype[ant]['genescores'] = {'wzx': otype[oantigen]["genescores"]["wzx"]}
            elif 'wzy' in otype[oantigen]["genescores"].keys():
                serotype[ant]['genescores'] = {'wzy': otype[oantigen]["genescores"]["wzy"]}
            selectedOantigen = oantigen
            break

        elif 'wzm' in otype[oantigen]["genescores"].keys() or 'wzt' in otype[oantigen]["genescores"].keys():
            serotype['O']['serogroup'] = oantigen
            if 'wzm' in otype[oantigen]["genescores"].keys():
                serotype[ant]['genescores'] = {'wzm': otype[oantigen]["genescores"]["wzm"]}
            elif 'wzt' in otype[oantigen]["genescores"].keys():
                serotype[ant]['genescores'] = {'wzt': otype[oantigen]["genescores"]["wzt"]}
            selectedOantigen = oantigen
            break


    if selectedOantigen:
        for allele in otype[selectedOantigen]["alleles"]:
            serotype[ant]["alleles"][allele] = setAlleleMeta(allelefieldnames, blastresultsdict[ant][allele])
            gene = otype[selectedOantigen]["allele2gene"][allele]
            serotype["gene2allele"][gene] = allele

    return serotype

def setAlleleMeta(keys,metadf):
    metaalleledict = {}
    for key in keys:
        metaalleledict[key] = metadf[key]
    return metaalleledict

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
        antigens = ["O","H"]
        for antigen in antigens:
            alleles = ectyper_dict[antigen]
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

    Otype = "-"; Htype = "-"

    numalleles = 0
    if "O" in final_results_dict[sample]:
        Otype=final_results_dict[sample]["O"]
        numalleles += len(final_results_dict[sample]["H"]["alleles"].keys())
    if "H" in final_results_dict[sample]:
        Htype=final_results_dict[sample]["H"]
        numalleles += len(final_results_dict[sample]["O"]["alleles"].keys())

    species=final_results_dict[sample]["species"]
    QCflag = "FAIL"

    if not re.match("Escherichia coli", species):
        return {"QCflag":QCflag,"NumberOfAlleles":numalleles}

    if Otype == "-" or Htype == "-":
        QCflag = "FAIL"
    else:
        QCflag = "PASS"

    return {"QCflag":QCflag,"NumberOfAlleles":numalleles}


def report_result(final_dict, output_dir, output_file):
    """
    Outputs the results of the ectyper run to the output file, and to the log.

    :param final_dict: Final ectyper predictions dictionary {'Sample': {'O': 'O26', 'H': 'H11', 'H11': {'fliC': 1.0}, 'O26': {'wzx': 0.46}}}
    :param output_file: File whose contents will be added to the log
    :return: None
    """

    header = "Name\tSpecies\tO-type\tH-type\tSerotype\tQC\t" \
             "Evidence\tGeneScores\tAllelesDBKeys\tQueryIdentities(%)\t" \
             "QueryCoverages(%)\tQueryContigNames\tQueryAlleleRanges\t" \
             "QueryAlleleLengths\tEcoliSpecificMarkers\tWarnings\n"
    output = []
    LOG.info(header.strip())

    alleleseqs={}
    for sample in final_dict.keys():
        output_line = [sample] #name of a query sample/genome
        output_line.append(final_dict[sample]["species"]) #add species info

        Otype = "-"; Htype = "-"
        if "O" in final_dict[sample].keys():
            Otype=final_dict[sample]["O"]["serogroup"]
        output_line.append(Otype)

        if "H" in final_dict[sample].keys():
            Htype=final_dict[sample]["H"]["serogroup"]
        output_line.append(Htype)

        output_line.append("{}:{}".format(Otype,Htype)) #serotype

        genescoresstr = ""; genelistorder = []
        if Otype != "-":
            for Ogene, score in sorted(final_dict[sample]["O"]["genescores"].items()): #makes sure that wzx/wzy and wzm/wzt are reported alphabetically
                genescoresstr = genescoresstr + "{}:{:.3f};".format(Ogene, score)
                genelistorder.append(("O", Ogene))

        if Htype != "-":
            for Hgene, score in final_dict[sample]["H"]['genescores'].items():
                genescoresstr = genescoresstr + "{}:{:.3f};".format(Hgene, score)
                genelistorder.append(("H", Hgene))


        QCdict = quality_control_results(sample, final_dict)
        output_line.append(QCdict["QCflag"])

        if genescoresstr == "":
            output_line.append("-")
        else:
            output_line.append("Based on {} allele(s)".format(QCdict["NumberOfAlleles"]))

        if genescoresstr == "":
            output_line.append("-")
        else:
            output_line.append(genescoresstr)

        identity_list = []; coverage_list = []; allele_list = []; contig_list = []
        startend_positions_list = []; length_match = []

        for ant, gene in genelistorder:
            allele = final_dict[sample]["gene2allele"][gene]
            identity_list.append(final_dict[sample][ant]["alleles"][allele]["identity"])
            coverage_list.append(final_dict[sample][ant]["alleles"][allele]["coverage"])
            contig_list.append(final_dict[sample][ant]["alleles"][allele]["contigname"])
            startend_positions_list.append(final_dict[sample][ant]["alleles"][allele]["startpos"]+"-"+
                                           final_dict[sample][ant]["alleles"][allele]["endpos"])
            length_match.append(final_dict[sample][ant]["alleles"][allele]["length"])
            allele_list.append(allele)

        for list in [allele_list,identity_list,coverage_list,contig_list,startend_positions_list,length_match]:
            if list:
                output_line.append(";".join(list))
            else:
                output_line.append("-")


        if "ecolimarkers" in final_dict[sample]:
            output_line.append(str(final_dict[sample]['ecolimarkers']))
        else:
            output_line.append("-")

        output_line.append(final_dict[sample]["error"])

        print_line = "\t".join(output_line)
        output.append(print_line + "\n")
        LOG.info(print_line)

        if alleleseqs:
            with open(file=output_dir+"/queryalleleseqs.fasta", mode="w") as fp:
                for allele, seq in alleleseqs.items():
                    fp.write(">{}\n{}\n".format(allele[1:],seq))
                fp.close()

    with open(output_file, mode="w", encoding="utf8") as ofh:
        ofh.write(header)
        for line in sorted(output):
            ofh.write(line)


def add_non_predicted(all_genomes_list, predictions_dict, other_dict, filesnotfound_dict, ecoli_dict):
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
            if gname in other_dict:
                predictions_dict[gname] = {
                    'error': other_dict[gname]["error"],
                    'species': other_dict[gname]["species"]
                }
            elif gname in filesnotfound_dict:
                predictions_dict[gname] = {
                'error': filesnotfound_dict[gname]["error"],
                'species': '-'
            }
            else:
                predictions_dict[gname] = {
                    'error':  "No O and H antigen determinant E.coli genes were found. Try running with --verify parameter",
                    'species': ecoli_dict[gname]["species"]
                }

    return predictions_dict
