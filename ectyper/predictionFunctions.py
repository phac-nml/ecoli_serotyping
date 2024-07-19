#!/usr/bin/env python

import json
import logging
import os, re
import pandas as pd
import ectyper.definitions as definitions
from ectyper import subprocess_util

from Bio import SeqRecord
from Bio.Seq import  Seq
from Bio.SeqIO import FastaIO

LOG = logging.getLogger(__name__)

"""
    Serotype prediction for E. coli
"""
def fasta_custom_header_generator(record):
    """Generate custom faster header for a fasta file using ID|Name|Description supplied to function

    Args:
        record (Bio.SeqRecord): BioPython Sequence object

    Returns:
        string: formatted FASTA header string using the format <id>|<accession>|<marker>|<gene symbol>|<gene full name>
    """
    return f"{record.id}|{record.name}|{record.description}"

def load_json(json_file):
    with open(json_file) as fp:
        return json.load(fp)

def json2fasta(json_file, output_dir):
    """

    """
    basename = os.path.basename(definitions.PATHOTYPE_ALLELE_JSON).split('.json')[0]
    fasta_pathotypedb_path = os.path.join(output_dir,f'{basename}.fasta')
    
    if os.path.exists(fasta_pathotypedb_path):
        return fasta_pathotypedb_path
    
    with open(json_file) as fp:
        json_db = json.load(fp)
    sequences = []
    
    for marker in json_db['markers'].keys():
        marker_dict = json_db['markers'][marker]
        for item in marker_dict:
            seq_name = item['accession']+'|'+marker+'|'+item['gene']
            if 'subtype' in item and 'stx' in marker:
                seq_name = seq_name+"|"+item['subtype']
            sequences.append(SeqRecord.SeqRecord(seq= Seq(item['dnasequence']) , id = item['id'], name= seq_name, 
                                                 description = re.sub(' ','_',item["description"]))) 
    with open(fasta_pathotypedb_path, "w") as fp:
        fasta_out = FastaIO.FastaWriter(fp, wrap=None, record2title=fasta_custom_header_generator )
        fasta_out.write_file(sequences)
    LOG.info(f"Created pathotype database in fasta format from JSON at {fasta_pathotypedb_path}") 
    return fasta_pathotypedb_path

def shiga_toxing_subtyping(pathotype_genes_tmp_df, output_dir, debug):
    results_dict = {
                    'stx_genes':[],
                    'stx_accessions':[],
                    'stx_allele_ids':[],
                    'stx_pidents': [],
                    'stx_pcovs':[],
                    'stx_contigs':[],
                    'stx_contigs_num': [],
                    'stx_gene_lengths':[],
                    'stx_gene_ranges':[]
                }
    for gene in ['stx1', 'stx2']:
        stx_toxin_df = pathotype_genes_tmp_df.query(f'gene == "{gene}"')
        stx_toxin_df.loc[:,['stx_subtype', 'stx_seqtype', 'stx_contig_name','rangeid']] = '-'
   
          
        
        if stx_toxin_df.empty == False:   
            #homogenize allele coordinates for +1 and -1 strands making coordinates on +1 strand
            hits_coordinates = stx_toxin_df.query('sframe == -1')[['send','sstart']]
            stx_toxin_df.loc[list(stx_toxin_df.loc[:,'sframe'] == -1),'sstart'] = hits_coordinates['send']
            stx_toxin_df.loc[list(stx_toxin_df.loc[:,'sframe'] == -1), 'send'] = hits_coordinates['sstart']
            #annotate stx dataframe of potential hits with shiga toxin subtype, contig name and sequence type
            stx_toxin_df.loc[:,'stx_subtype'] = stx_toxin_df['qseqid'].apply(lambda x:x.split('|')[4]).copy().to_list() #extract subtypes
            stx_toxin_df.loc[:,'stx_contig_name']  = stx_toxin_df['sseqid'].apply(lambda x: x.split('|')[2]).copy().to_list() #extract contig names
            stx_contigs = stx_toxin_df.stx_contig_name.unique()
            results_dict['stx_contigs_num'].append(f"{gene}:{len(stx_contigs)}")
            stx_types = ['complete_sequence','subunit_A', 'subunit_B'] #extract molecule types
            stx_seqtype = stx_toxin_df['qseqid'].apply(lambda x: stx_types[[i for i, type in enumerate(stx_types) if type in x][0]]).to_list()
            stx_toxin_df.loc[:,'stx_seqtype'] = stx_seqtype 

            # for each contig cluster stx ranges by based on coordinates overlap. create groups of overlapping ranges
            clusters_of_ranges = {}; cluster_counter=0
            for contig_name in stx_contigs:
                stx_toxin_df_tmp = stx_toxin_df.query(f'stx_contig_name == "{contig_name}"')
                unique_ranges = sorted(list(stx_toxin_df_tmp.apply(lambda x: range(x['sstart'],x['send']+1), axis=1).unique()),key=lambda r: r.start)
                ranges_groupby_overlap=[]
                #for each range find with which ranges it might overlap (all vs all)
                for urange in unique_ranges:
                    overlap = [True if len(set(urange).intersection(unique_ranges[i])) > 0 else False for i in range(0, len(unique_ranges))]
                    ranges_groupby_overlap.append([unique_ranges[i] for i,v in enumerate(unique_ranges) if overlap[i] == True])
                ranges_grouped_list = list(set(tuple(i) for i in ranges_groupby_overlap))
                
                #for each group of overlapping ranges create single range based on min and max values of ranges in that group
                ranges_candidate_list = []
                for range_list in ranges_grouped_list:
                    range_min = min([i[0] for i in range_list ]); range_max = max([i[-1] for i in range_list])+1
                    ranges_candidate_list.append(range(range_min,range_max))
                ranges_candidate_set_unique = [(i,len(i)) for i in set(ranges_candidate_list)]
                #final dictionary result of unique non-overlapping ranges per contig stored in clusters_of_ranges
                for range_index, (range_unique, range_length ) in enumerate(ranges_candidate_set_unique):
                    #find if a given range still overlaps with any in the set (exclude overlap of ranges with itself)
                    overlap_bool = [True if len(set(r).intersection(range_unique)) > 0 and (range_unique,l) not in ranges_candidate_set_unique else False for r, l in ranges_candidate_set_unique]
                    if all(overlap_bool) == False: #if no overlap is found select that range as final or else continue selection on max range length
                        clusters_of_ranges[f"{cluster_counter}"]={'range': range_unique, 'contig':contig_name, 'range_length': range_length}
                    else:
                        overlap_indices = [i for i, b in enumerate(overlap_bool) if b == True and i != range_index] 
                        overlap_ranges = [r for i,r in enumerate(ranges_candidate_set_unique) if i in overlap_indices]
                        overlap_max_len_range = max([r[1] for i,r in enumerate(overlap_ranges) if i in overlap_indices])
                        final_select_range = [r for r, l in ranges_candidate_set_unique if l == overlap_max_len_range][0]
                        clusters_of_ranges[f"{cluster_counter}"]={'range': final_select_range, 'contig':contig_name, 'range_length': range_length}
                    cluster_counter += 1       
            #annotate shiga toxin dataframe with ranges identifiers
            stx_toxin_df.loc[:,'rangeid'] = stx_toxin_df.apply(lambda x: [range_id 
                                                                        for range_id in clusters_of_ranges.keys() 
                                                                        if x['sstart'] in clusters_of_ranges[range_id]['range'] and
                                                                        x['stx_contig_name'] == clusters_of_ranges[range_id]['contig']][0], 
                                                                        axis=1).copy().to_list()
            stx_toxin_df = stx_toxin_df.astype({'rangeid':'int32', 'bitscore':'int32'})
            LOG.info(f"Total of {len(clusters_of_ranges)} non-overlapping ranges found on {len([v['contig'] for v in clusters_of_ranges.values()])} contigs for {gene} with ranges ({[ {v['contig']: (min(v['range']),max(v['range']))} for v in clusters_of_ranges.values()]} )".rstrip("\n"))
            #write shiga toxin dataframe for inspection just in case
            if debug:
                stx_df_out_filename = f'{gene}_allhits_annotated_df.txt'
                LOG.debug(f"Wrote {gene} annotated potential hits dataframe to {output_dir}/{stx_df_out_filename}")
                stx_toxin_df.to_csv(os.path.join(output_dir,stx_df_out_filename), sep="\t", index=False)
            # get top hit for each common gene range. Provide mixed call if >1 hits share the same 'bitscore'
            stx_subtypes_dict={}
            for range_id in stx_toxin_df['rangeid'].unique():
                max_bitscore = stx_toxin_df.query(f'rangeid == {range_id}')['bitscore'].max()
                tmp_df = stx_toxin_df.query(f'rangeid == {range_id} and bitscore == {max_bitscore}')
                for subtype in tmp_df['stx_subtype'].unique():
                    stx_hit = tmp_df.query(f'`stx_subtype` == "{subtype}"').iloc[0]
                    if subtype not in stx_subtypes_dict: #insert only subtypes not already encountered in the dictionary even if that would a data loss as usually those a lower quality results
                        stx_subtypes_dict[subtype] = stx_hit[['sstart','send']].to_dict()
                        stx_subtypes_dict[subtype]['index'] = stx_hit.name
                        stx_subtypes_dict[subtype]['bitscore'] = stx_hit.bitscore
                        stx_subtypes_dict[subtype]['contig'] = stx_hit.stx_contig_name
                        stx_subtypes_dict[subtype]['length'] = stx_hit.length
                        stx_subtypes_dict[subtype]['gene_range'] = f"{stx_hit.sstart}-{stx_hit.send}"
                        stx_subtypes_dict[subtype]['gene'] = stx_hit.gene+stx_hit.stx_subtype
                        stx_subtypes_dict[subtype]['accession'] = stx_hit.accession
                        stx_subtypes_dict[subtype]['allele_id'] = stx_hit.allele_id
                        stx_subtypes_dict[subtype]['pident'] = stx_hit.pident
                        stx_subtypes_dict[subtype]['pcov'] = stx_hit.qcovhsp
                        stx_subtypes_dict[subtype]['subtype'] = subtype
                       
        
            for contig in set([stx_subtypes_dict[s]['contig'] for s in stx_subtypes_dict.keys()]):
                for dict_item in [stx_subtypes_dict[s] for s in stx_subtypes_dict if stx_subtypes_dict[s]['contig'] == contig]:
                    results_dict['stx_genes'].append(dict_item['gene'])
                    results_dict['stx_accessions'].append(dict_item['accession'])
                    results_dict['stx_allele_ids'].append(dict_item['allele_id'])
                    results_dict['stx_pidents'].append(str(dict_item['pident']))
                    results_dict['stx_pcovs'].append(str(dict_item['pcov']))
                    results_dict['stx_gene_lengths'].append(str(dict_item['length']))
                    results_dict['stx_contigs'].append(str(dict_item['contig']))
                    results_dict['stx_gene_ranges'].append(str(dict_item['gene_range']))
                
    #report shiga toxin subtypes in alphabetical order    
    sorted_order = [i[0] for i in sorted(enumerate(results_dict['stx_genes']), key=lambda x:x[1])  ] 
    for k in results_dict.keys():
        if len(results_dict[k]) != 0:
            results_dict[k] = ";".join([results_dict[k][i] for i in sorted_order])
        else:
            results_dict[k] = "-"   
    print(results_dict)        
    return results_dict

def predict_pathotype_and_shiga_toxin_subtype(ecoli_genome_files_dict, other_genomes_dict, temp_dir, verify_species_flag, pident, pcov, 
                                              output_dir, debug, pathotype_db):
    """Get pathotype(s) of a sample

    Args:
        ecoli_genome_files_dict (dict): dictionary of E.coli genome paths to analyze (key: sample_id, value: path)
        other_genomes_dict (dict): dictionary of non-E.coli genome paths to analyze (key: sample_id, value: path)
        temp_dir (path): temporary folder path
        pident (float): min %identity to filter blastn pathotype allele hits
        pcov (float): min %coverage to filter blastn pathotype allele hits
        output_dir (path): output results directory path

    Returns:
        list: list of pathotypes
    """
    LOG.info(f"Starting pathotype predictions on {len(ecoli_genome_files_dict.keys())} samples. Reminder: Please use --verify option to run pathotype predictions only on E.coli samples ...")

    if len(other_genomes_dict.keys()) > 0 and verify_species_flag == True:
        LOG.info(f"A total of {len(other_genomes_dict.keys())} non-E.coli sample(s) will not be pathotyped. Omit --verify option if you want to type ALL samples regardless ...")
    
    path2patho_db = json2fasta(definitions.PATHOTYPE_ALLELE_JSON, temp_dir)
    json_patho_db = load_json(definitions.PATHOTYPE_ALLELE_JSON)
    
    predictions_pathotype_dict = {}
    pathotype_genes_overall_df = pd.DataFrame()

    #perform pathotyping for all samples if --pathotype key is provided unless a --verify option is specified
    merged_samples_dict = ecoli_genome_files_dict.copy()
    merged_samples_dict.update(other_genomes_dict)

    for g in merged_samples_dict.keys():
        LOG.info(f"Running pathotype prediction on {g} ...")
        predictions_pathotype_dict[g]={field:'-' for field in definitions.PATHOTYPE_TOXIN_FIELDS}
        if verify_species_flag == True and 'Escherichia coli' not in merged_samples_dict[g]['species']:
            LOG.info(f"Skipping pathotype prediction for {g} as species is not E.coli and --verify is specified ...")
            continue
        input_sequence_file = merged_samples_dict[g]['modheaderfile'] 
        cmd = [
            "blastn",
            "-query", path2patho_db,
            "-subject", input_sequence_file,
            "-out", f"{temp_dir}/blast_pathotype_result.txt",
            "-outfmt", "6 qseqid qlen sseqid length pident sstart send sframe slen qcovhsp bitscore sseq"
        ]
        LOG.debug(f"BLASTN results on pathotype database written to {temp_dir}/blast_pathotype_result.txt ...")
        subprocess_util.run_subprocess(cmd)
        if os.stat(f'{temp_dir}/blast_pathotype_result.txt').st_size == 0:
            LOG.warning(f"No pathotype signatures found for sample {g} as pathotype BLAST results file {temp_dir}/blast_pathotype_result.txt is empty. Skipping ...")
            predictions_pathotype_dict[g]={field:'-' for field in definitions.PATHOTYPE_TOXIN_FIELDS}
            predictions_pathotype_dict[g]['pathotype']= ['ND']
            continue
        
        pathotype_genes_tmp_df = pd.read_csv(f'{temp_dir}/blast_pathotype_result.txt', sep="\t", header=None)
        pathotype_genes_tmp_df.columns = ['qseqid', 'qlen', 'sseqid', 'length', 'pident', 'sstart', 'send', 'sframe', 'slen', 'qcovhsp', 'bitscore', 'sseq']
        pathotype_genes_tmp_df.sort_values('bitscore', ascending=False,inplace=True) #sort hit in descending order by bitscore
        pathotype_genes_tmp_df['allele_id'] = pathotype_genes_tmp_df['qseqid'].apply(lambda x:x.split('|')[0])
        pathotype_genes_tmp_df['accession'] = pathotype_genes_tmp_df['qseqid'].apply(lambda x:x.split('|')[1])
        pathotype_genes_tmp_df['gene'] = pathotype_genes_tmp_df['qseqid'].apply(lambda x:x.split('|')[2])
        pathotype_genes_tmp_df['sample_id'] = g
        
        pathotype_genes_tmp_df.query(f'pident >= {pident} and qcovhsp >= {pcov}', inplace=True) # Default BioNumerics threshold in pathotype module

        #shiga toxin subtyping
        result_stx1_stx2 = shiga_toxing_subtyping(pathotype_genes_tmp_df, output_dir, debug)
        predictions_pathotype_dict[g].update(result_stx1_stx2) #Update dictionary with shiga toxin typing results
        
        
        if debug == True:
            pathotype_genes_overall_df = pd.concat([pathotype_genes_overall_df,pathotype_genes_tmp_df], ignore_index=True)
        pathotype_genes_top_hits = pathotype_genes_tmp_df.loc[pathotype_genes_tmp_df.groupby('gene')['bitscore'].transform("idxmax").unique()].sort_values('gene')
        pathotype_genes_top_hits = pathotype_genes_top_hits.sort_values('gene', axis=0)
        
        gene_list = pathotype_genes_top_hits['gene'].to_list()
        accessions_list = pathotype_genes_top_hits['accession'].to_list()
        if len(gene_list)!= 0:
            predictions_pathotype_dict[g]['pathotype'] = []
            predictions_pathotype_dict[g]['pathotype_rule_ids']=[]
            predictions_pathotype_dict[g]['pathotype_genes'] = ",".join(gene_list)
            predictions_pathotype_dict[g]['pathotype_gene_counts'] = {}
            predictions_pathotype_dict[g]['pathotype_gene_names'] = ",".join([f"{g}: {r['description']}" for g, a in zip(gene_list,accessions_list) for r  in json_patho_db['markers'][g] if r['accession'] == a])
            predictions_pathotype_dict[g]['pathotype_allele_id']=";".join(pathotype_genes_top_hits['allele_id'].to_list())
            predictions_pathotype_dict[g]['pathotype_accessions']=";".join(accessions_list)
            predictions_pathotype_dict[g]['pathotype_pident'] = ";".join(pathotype_genes_top_hits['pident'].astype(str).to_list())
            predictions_pathotype_dict[g]['pathotype_pcov'] = ";".join(pathotype_genes_top_hits['qcovhsp'].astype(str).to_list())
            predictions_pathotype_dict[g]['pathotype_length_ratio'] = ";".join(pathotype_genes_top_hits.apply(lambda row : f"{row['length']}/{row['qlen']}", axis=1).to_list())
            #pathotypes assigner logic
            for pathotype in pathotype_db['pathotypes']:
                pathotype_genes_found = []
                for rule_id in pathotype_db['pathotypes'][pathotype]['rules'].keys():
                    ngenes_in_rule = len(pathotype_db['pathotypes'][pathotype]['rules'][rule_id]['genes'])
                    matched_gene_counter = 0
                    for gene in pathotype_db['pathotypes'][pathotype]['rules'][rule_id]['genes']:
                        if "!" in gene:
                            gene = re.sub('!','',gene)
                            if gene not in gene_list:
                                matched_gene_counter += 1
                            else:
                                matched_gene_counter = 0
                                break    
                        elif gene in gene_list:
                            pathotype_genes_found.append(gene)
                            matched_gene_counter += 1
                    #if pathotype rule was completely matched        
                    if ngenes_in_rule == matched_gene_counter:
                        predictions_pathotype_dict[g]['pathotype'].append(pathotype)
                        predictions_pathotype_dict[g]['pathotype_rule_ids'].extend([f"{rule_id}:{p}:{'+'.join(json_patho_db['pathotypes'][p]['rules'][r]['genes'] )}" for p  in json_patho_db['pathotypes'] for r in json_patho_db['pathotypes'][p]['rules'] if r == rule_id])              
                        predictions_pathotype_dict[g]['pathotype_gene_counts'][pathotype] = []
                        predictions_pathotype_dict[g]['pathotype_gene_counts'][pathotype].extend(pathotype_genes_found) 

        final_pathotypes_list = list(set(predictions_pathotype_dict[g]['pathotype']))
        if 'EHEC' in final_pathotypes_list:
            final_pathotypes_list.remove('EHEC'); final_pathotypes_list.remove('STEC')
            final_pathotypes_list.append('EHEC-STEC')
        if '-' in final_pathotypes_list  or len(final_pathotypes_list) == 0:
            predictions_pathotype_dict[g]['pathotype'] = ['ND']
            predictions_pathotype_dict[g]['pathotype_count'] = "0"
            predictions_pathotype_dict[g]['pathotype_rule_ids'] = '-'
            predictions_pathotype_dict[g]['pathotype_gene_counts']= '-'
        else:    
            predictions_pathotype_dict[g]['pathotype'] = final_pathotypes_list #final pathotype(s) list
            predictions_pathotype_dict[g]['pathotype_count'] = f"{len(final_pathotypes_list)}"
            predictions_pathotype_dict[g]['pathotype_rule_ids'] = ";".join(predictions_pathotype_dict[g]['pathotype_rule_ids'])
            predictions_pathotype_dict[g]['pathotype_gene_counts']= ";".join(sorted([f"{p}:{len(set(predictions_pathotype_dict[g]['pathotype_gene_counts'][p]))} ({','.join(sorted(set(predictions_pathotype_dict[g]['pathotype_gene_counts'][p])))})" for p in predictions_pathotype_dict[g]['pathotype_gene_counts']]))
        # add pathotype database version
        predictions_pathotype_dict[g]['pathotype_database'] = f"v{json_patho_db['version']}" 
    #write pathotype blastn results
    if debug == True:
        LOG.debug(f"Writting overall pathotype BLASTn results to {output_dir}/blastn_pathotype_alleles_overall.txt")
        pathotype_genes_overall_df.to_csv(f'{output_dir}/blastn_pathotype_alleles_overall.txt',sep="\t", index=False)
    

    return predictions_pathotype_dict
        

def predict_serotype(blast_output_file, ectyper_dict, args):
    """
    Predict the serotype of all genomes, given the blast output of the markers against the genomes

    :param blast_output_file: BLAST results of O and H-type allele search against the genomes of interest
    :param ectyper_dict: ectyper database in dict format from JSON file of known alleles and their O and H mappings
    :param args: Commandline arguments
    :return: The CSV formatted predictions file
    """

    LOG.info(f"Predicting serotype from BLAST output {blast_output_file}")
    if os.stat(blast_output_file).st_size == 0:
        LOG.warning(f"Empty O- and H- antigen BLAST output file {blast_output_file}. O-/H- antigen prediction will be skipped as no potential hits available. Check inputs")
        return {}, pd.DataFrame()
        
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
        predictions_dict[genome_name] = get_prediction(per_genome_df)
    LOG.info("Serotype prediction successfully completed")
    LOG.debug("Predictions dict:\n{}".format(predictions_dict))

    return predictions_dict, output_df


def get_prediction(per_genome_df):
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
            blastresultsdict[ant][row.qseqid]["contigname"] = row.sseqid.split('|')[2] #3rd position is used for contig name assignment
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
    best_order_list = sorted([item[0] for item in scorestupleslist], key=lambda x: int(x[1:])) #['O102']


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

    
    header = "\t".join(definitions.OUTPUT_TSV_HEADER)+"\n"
    
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
        if 'pathotype' in final_dict[sample]:
            for field in definitions.PATHOTYPE_TOXIN_FIELDS:
                output_line.append(final_dict[sample][field])
        else:
            output_line.extend(['-']*len([f for f in definitions.OUTPUT_TSV_HEADER if 'Pathotype' in f or 'Stx' in f]))
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
