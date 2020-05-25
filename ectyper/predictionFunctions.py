#!/usr/bin/env python

import json
import logging
import os, re
import pandas as pd
import ectyper.definitions as definitions


LOG = logging.getLogger(__name__)

"""
    Serotype prediction for E. coli
"""


def predict_serotype(blast_output_file, ectyper_dict, args):
    """
    Predict the serotype of all genomes, given the blast output of the markers against the genomes

    :param blast_output_file: BLAST results of O and H-type allele search against the genomes of interest
    :param ectyper_dict: ectyper database in dict format from JSON file of known alleles and their O and H mappings
    :param args: Commandline arguments
    :return: The CSV formatted predictions file
    """

    LOG.info("Predicting serotype from blast output")
    output_df = blast_output_to_df(blast_output_file) #columns: length	pident	qcovhsp	qlen qseqid	send sframe	sseq	sseqid	sstart	score

    ectyper_df = ectyper_dict_to_df(ectyper_dict) #columns: antigen	desc gene name


    # Merge output_df and ectyper_df
    output_df = output_df.merge(ectyper_df, left_on='qseqid', right_on='name', how='left')
    predictions_dict = {}


    OantigensPotential = set(output_df.query('type == "O" & pident > '+str(definitions.MIN_O_IDENTITY_LS)
                                            +' & qcovhsp > '+str(definitions.MIN_O_COVERAGE_LS ))["antigen"])

    if len(OantigensPotential) == 1:
        output_df = output_df.query(
            '(type == "O" & pident >= ' + str(definitions.MIN_O_IDENTITY_LS) + ' & qcovhsp >= ' +
            str(definitions.MIN_O_COVERAGE_LS) + ') | '
            '(type == "H" & pident >= ' +str(args.percentIdentityHtype) + ' & qcovhsp >= ' + str(args.percentCoverageHtype) + ' )')
    else:
        output_df = output_df.query('(type == "O" & pident >= '+str(args.percentIdentityOtype)+' & qcovhsp >= '+str(args.percentCoverageOtype)+') | '
                                    '(type == "H" & pident >= '+str(args.percentIdentityHtype)+' & qcovhsp >= '+str(args.percentCoverageHtype)+' )')

    # Append individual genomes/samples to overall dataframe for future selection
    output_df = output_df.assign(genome_name=output_df['sseqid'].str.split('|').str[1].values)

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
    #serotype["gene2allele"] = {}


    # Go for the highest match, if both genes exist over the thresholds
    blastresultsdict={}; blastresultsdict["O"]={}; blastresultsdict["H"]={}

    for row in per_genome_df.itertuples():
        # H is already set, skip
        # get the 'O' or 'H' from the antigen column
        ant = row.antigen[:1]
        if row.qseqid not in blastresultsdict[ant]:  #row.qseqid = allele database key (e.g. H29-2-fliC-origin)
            blastresultsdict[ant][row.qseqid] = {}
            blastresultsdict[ant][row.qseqid]["gene"] = row.gene
            blastresultsdict[ant][row.qseqid]["antigen"] = row.antigen
            blastresultsdict[ant][row.qseqid]["score"] = float(row.score)
            blastresultsdict[ant][row.qseqid]["identity"] = float(row.pident)
            blastresultsdict[ant][row.qseqid]["coverage"] = float(row.qcovhsp)
            blastresultsdict[ant][row.qseqid]["contigname"] = row.sseqid.split('|')[2]
            blastresultsdict[ant][row.qseqid]["startpos"] = int(row.sstart)
            blastresultsdict[ant][row.qseqid]["endpos"] = int(row.send)
            blastresultsdict[ant][row.qseqid]["length"] = int(row.length)
            blastresultsdict[ant][row.qseqid]["shared"] = bool(row.sharedallele)

    allelefieldnames = ["identity", "coverage", "contigname", "length", "startpos", "endpos", "gene"]

    ant = "H"
    if blastresultsdict[ant].keys():
        topHallele = sorted([(dballele,blastresultsdict[ant][dballele]["score"]) for dballele in blastresultsdict["H"].keys()],
                            key=lambda x:x[1], reverse=True)[0][0]
        serotype[ant]["serogroup"] = blastresultsdict[ant][topHallele]["antigen"]
        gene = blastresultsdict[ant][topHallele]["gene"]
        score = blastresultsdict[ant][topHallele]["score"]
        serotype[ant]["genescores"] = {gene:score}
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



    # rank O-type serovars based on the sum of scores from BOTH alleles (wzx/wzy or wzt/wzm)
    # calculate score  for O antigen
    rank_Otype_dict = {}
    for oantigen in otype.keys():
        rank_Otype_dict[oantigen]={"numalleles": 0, "scores":[], "sumscore": 0}
        rank_Otype_dict[oantigen]["numalleles"] = len(otype[oantigen]["alleles"])
        if rank_Otype_dict[oantigen]["numalleles"] == 4:
            LOG.warning("O-antigen {} has 4 alleles instead of expected 2."
                        "Might be due extra genes pairs (wzx/wzy,wzt/wzm) mapping to the same antigen or multiple gene copies".format(oantigen))
            wzx_wzy_score_sum = otype[oantigen]["genescores"]['wzx'] + otype[oantigen]["genescores"]['wzy']
            wzm_wzt_score_sum = otype[oantigen]["genescores"]['wzm'] + otype[oantigen]["genescores"]['wzt']
            if wzx_wzy_score_sum > wzm_wzt_score_sum:
                rank_Otype_dict[oantigen]["scores"] = [otype[oantigen]["genescores"][gene] for gene in ['wzx','wzy']]
            else:
                rank_Otype_dict[oantigen]["scores"] = [otype[oantigen]["genescores"][gene] for gene in ['wzm', 'wzt']]

        elif rank_Otype_dict[oantigen]["numalleles"] == 3:
            if 'wzx' in otype[oantigen]["genescores"].keys() and 'wzy' in otype[oantigen]["genescores"].keys():
                rank_Otype_dict[oantigen]["scores"] = [otype[oantigen]["genescores"][gene] for gene in ['wzx', 'wzy']]
            elif 'wzm' in otype[oantigen]["genescores"].keys() and 'wzt' in otype[oantigen]["genescores"].keys():
                rank_Otype_dict[oantigen]["scores"] = [otype[oantigen]["genescores"][gene] for gene in ['wzm', 'wzt']]
        else:
            rank_Otype_dict[oantigen]["scores"] = otype[oantigen]["genescores"].values()
        rank_Otype_dict[oantigen]["sumscore"] = sum(rank_Otype_dict[oantigen]["scores"])



    scorestupleslist = [(otypename,rank_Otype_dict[otypename]["sumscore"]) for otypename in rank_Otype_dict.keys()]
    scorestupleslist = sorted(scorestupleslist, key=lambda x: x[1], reverse=True) #[('O102', 1.73)]
    best_order_list = [item[0] for item in scorestupleslist] #['O102']


    LOG.debug("Otype dict:{}".format(otype))
    LOG.debug("Serotype dict:{}".format(serotype))
    LOG.debug("Best order alleles-scores list of tuples:{}".format(scorestupleslist))
    LOG.debug("Best order alleles list:{}".format(best_order_list))

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


    #find pairwise distance between all alleles indentified by BLAST search and find alleles with IDENTICAL O-antigen score
    #Note that difference in score of 0 between O antigens is considered to be identical

    identicalscorestupleslist = [(orow, ocol, abs(i - j)) for ocol, i in scorestupleslist for orow, j in scorestupleslist
                            if i - j == 0 and (orow == selectedOantigen or ocol == selectedOantigen) and
                            (orow != ocol)]

    #find highly similar antigens at 99.9% indentity which is only 1 nt apart
    highsimilarity_oantigens = [orow for ocol, i in scorestupleslist for orow, j in
                                scorestupleslist if i - j <= definitions.HIGH_SIMILARITY_THRESHOLD_O
                                and orow != selectedOantigen
                                and (orow == selectedOantigen or ocol == selectedOantigen)
                                and (orow != ocol)]


    if highsimilarity_oantigens:
        mixedoantigen = [selectedOantigen] + highsimilarity_oantigens
        serotype['O']['serogroup'] = "/".join(mixedoantigen)
        LOG.info("Highly similar O-antigen candidates were found for {}".format(mixedoantigen))
    elif selectedOantigen != "":
        serotype['O']['serogroup'] = selectedOantigen

    # append to existing top O-antigen and generate mixed antigen call
    # alleles with identical score
    if identicalscorestupleslist:
        identical_oantigens = [i for i, j, s in identicalscorestupleslist if i != selectedOantigen]

        for oantigen in identical_oantigens:
            for allele in otype[oantigen]["alleles"]:
                serotype[ant]["alleles"][allele] = setAlleleMeta(allelefieldnames, blastresultsdict[ant][allele])


    # add info for selected alleles
    if selectedOantigen:
        for allele in otype[selectedOantigen]["alleles"]:
            serotype[ant]["alleles"][allele] = setAlleleMeta(allelefieldnames, blastresultsdict[ant][allele])


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
                'qlen':   fields[1],
                'sseqid': fields[2],
                'length': fields[3],
                'pident': fields[4],
                'sstart': fields[5],
                'send':   fields[6],
                'sframe': fields[7],
                'qcovhsp': fields[8],
                'bitscore': fields[9],
                'sseq': fields[10]
            }
            output_data.append(entry)
    df = pd.DataFrame(output_data)

    typeconvertdict={'qseqid':str,'qlen':int,'sseqid':str,'length':int,'pident':float,'sstart':int,'send':int,
                     'sframe':int, 'qcovhsp':float,'bitscore':float, 'sseq':str}
    df =  df.astype(typeconvertdict)

    if not output_data:
        LOG.warning("No hits found for blast output file {}".format(blast_output_file))

        # Return empty dataframe with correct columns
        return pd.DataFrame(
            columns=[
                'length', 'pident', 'qcovhsp',
                'qlen', 'qseqid', 'send',
                'sframe', 'sseqid', 'sstart','bitscore', 'sseq'
            ])
    else:
        df['score'] = df['pident'].astype(float)*df['qcovhsp'].astype(float)/10000
        return df


def ectyper_dict_to_df(ectyper_dict):
    """
    Load all the known alleles for the O and H genes, and store them as a DataFrame.
    :param ectyper_dict: ECTyper database as a dictionary from a JSON file of all known O and H alleles
    :return: DataFrame of the JSON file
    """

    temp_list = []
    antigens = ["O", "H"]
    for antigen in antigens:
        alleles = ectyper_dict[antigen]

        for name, allele in alleles.items():
            new_entry = {
                'type': antigen,
                'antigen': allele.get('allele'),
                'name': name,
                'gene': allele.get('gene'),
                'desc': allele.get('desc'),
                'sharedallele': allele.get('isAlleleShared')
            }
            temp_list.append(new_entry)
    df = pd.DataFrame(temp_list)

    return df

def mean(numbers):
    return sum(numbers)/len(numbers)

def getPredictionNumAlleles(sample, final_results_dict):
    numalleles = 0
    if "O" in final_results_dict[sample]:
        numalleles += len(final_results_dict[sample]["O"]["genescores"].keys()) #some O-antigens like O8 have both wzx/wzy and wzm/wzt
    if "H" in final_results_dict[sample]:
        numalleles += len(final_results_dict[sample]["H"]["genescores"].keys())
    return numalleles


def getQuality_control_results(sample, final_results_dict, ectyperdb_dict):
    """
    Determined approximate quality of the prediction based on the allele scores. Adopt pessimistic approach by looking at min values
    :param sample: sample/genome name)
    :param final_results_dict: dictionary with final output results (e.g. serovar, sequences, conf. scores)
    :return: Quality dictionary
    """

    if 'O' in final_results_dict[sample]:
        Otype = final_results_dict[sample]['O']['serogroup']
        Otypealleles = final_results_dict[sample]['O']['alleles'].keys()
    else:
        Otype = "-"
        Otypealleles=[]

    if 'H' in final_results_dict[sample]:
        Htype = final_results_dict[sample]['H']['serogroup']
        Htypealleles = final_results_dict[sample]['H']['alleles'].keys()
    else:
        Htype = "-"
        Htypealleles=[]

    #QC1: Check if species is E.coli
    if not re.match("Escherichia coli", final_results_dict[sample]["species"]):
        final_results_dict[sample]["error"]=final_results_dict[sample]["error"]+"Sample was not identified as valid E.coli sample but as {}".format(final_results_dict[sample]["species"])
        return "WARNING (WRONG SPECIES)"

    if not 'O' and 'H' in final_results_dict[sample]:
        return  "-"



    #QC2: Check on quality of O and H antigen prediction if at all taking into account resolution power of each antigen
    if Otype == "-" and Htype == "-":
        final_results_dict[sample]["error"]=final_results_dict[sample]["error"]+"Failed to type both O and H antigens. Consider lowering thresholds or traditional serotyping"
        return "FAIL (-:- TYPING)"
    elif Otype == "-" and Htype != "-":
        final_results_dict[sample]["error"]=final_results_dict[sample]["error"]+"Failed to type O antigen. Consider lowering thresholds or traditional serotyping"
        return "WARNING (-:H TYPING)"
    elif Otype != "-" and Htype == "-":
        final_results_dict[sample]["error"] = final_results_dict[sample]["error"] + "Failed to type H antigen. Consider lowering thresholds or traditional serotyping"
        return "WARNING (O:- TYPING)"
    elif len(Otype.split("/")) >= 2:
        final_results_dict[sample]["error"] = final_results_dict[sample]["error"] + "Mixed O-antigen call reported due to very high degree of antigen groups identity in excess of {}%." \
                                                                                    "Consider traditional serotyping as in-silico predictions might not be accurate".format((1-definitions.HIGH_SIMILARITY_THRESHOLD_O)*100)
        return "WARNING MIXED O-TYPE"
    else:

        checkpredscoresboolO = [final_results_dict[sample]["O"]["alleles"][allele]["identity"] *
                                final_results_dict[sample]["O"]["alleles"][allele]["coverage"] / 1e4 >=
                                ectyperdb_dict["O"][allele]["MinPident"] * ectyperdb_dict["O"][allele]["MinPcov"] / 1e4
                                for allele in Otypealleles]

        checkpredscoresboolH = [final_results_dict[sample]["H"]["alleles"][allele]["identity"] *
                                final_results_dict[sample]["H"]["alleles"][allele]["coverage"] / 1e4 >=
                                ectyperdb_dict["H"][allele]["MinPident"] * ectyperdb_dict["H"][allele]["MinPcov"] / 1e4
                                for allele in Htypealleles]

        if all(checkpredscoresboolO+checkpredscoresboolH) and all([item != [] for item in [checkpredscoresboolO,checkpredscoresboolH]]):
            QCflag="PASS (REPORTABLE)"
        elif all(checkpredscoresboolO):
            QCflag="WARNING (H NON-REPORT)"
            final_results_dict[sample]["error"]=final_results_dict[sample]["error"]+"H-antigen has %identity or %coverage below reportable threshold"
        elif all(checkpredscoresboolH):
            QCflag="WARNING (O NON-REPORT)"
            final_results_dict[sample]["error"]=final_results_dict[sample]["error"]+"O-antigen has %identity or %coverage below reportable threshold"
        else:
            QCflag="WARNING (O and H NON-REPORT)"
            final_results_dict[sample]["error"] = final_results_dict[sample]["error"] + "O and H-antigens have %identity or %coverage below reportable threshold"


    return QCflag



def report_result(final_dict, output_dir, output_file, args):
    """
    Outputs the results of the ectyper run to the output file, and to the log.

    :param final_dict: Final ectyper predictions dictionary {'Sample': {'O': 'O26', 'H': 'H11', 'H11': {'fliC': 1.0}, 'O26': {'wzx': 0.46}}}
    :param output_file: File whose contents will be added to the log
    :return: None
    """

    header = "Name\tSpecies\tO-type\tH-type\tSerotype\tQC\t" \
             "Evidence\tGeneScores\tAlleleKeys\tGeneIdentities(%)\t" \
             "GeneCoverages(%)\tGeneContigNames\tGeneRanges\t" \
             "GeneLengths\tDatabase\tWarnings\n"
    output = []
    LOG.info(header.strip())

    alleleseqs={}
    for sample in final_dict.keys():
        output_line = [sample] #name of a query sample/genome
        output_line.append(final_dict[sample]["species"]) #add species info

        if "O" in final_dict[sample].keys():
            Otype=final_dict[sample]["O"]["serogroup"]
        else:
            Otype = "-"
        output_line.append(Otype)

        if "H" in final_dict[sample].keys():
            Htype=final_dict[sample]["H"]["serogroup"]
        else:
            Htype = "-"
        output_line.append(Htype)

        output_line.append("{}:{}".format(Otype,Htype)) #serotype

        genescoresstr = ""; geneordertuplelist = []; genelistordersimple =[]
        if Otype != "-":
            for Ogene, score in sorted(final_dict[sample]["O"]["genescores"].items()): #makes sure that wzx/wzy and wzm/wzt are reported alphabetically
                genescoresstr = genescoresstr + "{}:{:.3g};".format(Ogene, score)
                geneordertuplelist.append(("O", Ogene, score))
                genelistordersimple.append(Ogene)

        if Htype != "-":
            for Hgene, score in final_dict[sample]["H"]['genescores'].items():
                genescoresstr = genescoresstr + "{}:{:.2g};".format(Hgene, score)
                geneordertuplelist.append(("H", Hgene, score))
                genelistordersimple.append(Hgene)




        if args.verify is True and 'QC' in final_dict[sample]:
            QCflag = final_dict[sample]["QC"]
            output_line.append(QCflag)
        else:
            output_line.append("-")


        if genescoresstr == "":
            output_line.append("-")
        else:
            output_line.append("Based on {} allele(s)".format(getPredictionNumAlleles(sample,final_dict)))

        if genescoresstr == "":
            output_line.append("-")
        else:
            output_line.append(genescoresstr)


        reportdict={}
        for ant, gene, score in geneordertuplelist :
            reportdict[gene]={"identity":[],"coverage":[],"alleles":[],"contigname":[],"coordinates":[], "length":[]}
            for allele, alleledict in final_dict[sample][ant]["alleles"].items():
                if final_dict[sample][ant]["alleles"][allele]["gene"] == gene:
                    reportdict[gene]["identity"].append("{:g}".format(final_dict[sample][ant]["alleles"][allele]["identity"]))
                    reportdict[gene]["coverage"].append("{:g}".format(final_dict[sample][ant]["alleles"][allele]["coverage"]))
                    reportdict[gene]["alleles"].append(allele)
                    reportdict[gene]["contigname"].append(final_dict[sample][ant]["alleles"][allele]["contigname"])
                    reportdict[gene]["length"].append("{}".format(final_dict[sample][ant]["alleles"][allele]["length"]))
                    if final_dict[sample][ant]["alleles"][allele]["startpos"] < final_dict[sample][ant]["alleles"][allele]["endpos"]:
                        reportdict[gene]["coordinates"].append("{}-{}".format(final_dict[sample][ant]["alleles"][allele]["startpos"],
                                                                              final_dict[sample][ant]["alleles"][allele]["endpos"]))
                    else:
                        reportdict[gene]["coordinates"].append("{}-{}".format(final_dict[sample][ant]["alleles"][allele]["endpos"],
                                                                              final_dict[sample][ant]["alleles"][allele]["startpos"]))
            for field in ["identity","coverage","contigname","coordinates","length"]:
                reportdict[gene][field] = set(reportdict[gene][field])



        for field in ["alleles", "identity","coverage","contigname","coordinates","length"]:
            finalvalue = ""
            for ant, gene, score in geneordertuplelist:
                if field == "alleles" and ant=="O":
                    valuesforgene = [item for o in Otype.split("/") for item in reportdict[gene]["alleles"] if re.match(o, item)]
                else:
                    valuesforgene = reportdict[gene][field]

                if valuesforgene:
                    finalvalue = finalvalue+"/".join(valuesforgene)+";"
                else:
                    finalvalue = finalvalue+"-;"
            output_line.append(finalvalue)

        output_line.append(final_dict[sample]["database"])

        if final_dict[sample]["error"]:
            output_line.append(final_dict[sample]["error"])
        else:
            output_line.append("-")

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
