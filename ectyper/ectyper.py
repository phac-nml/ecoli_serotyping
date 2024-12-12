#!/usr/bin/env python
"""
    Predictive serotyping for _E. coli_.
"""
import os, shutil, random
import datetime
import json, gzip
import logging
from multiprocessing import Pool
from functools import partial

from ectyper import (commandLineOptions, definitions, speciesIdentification,
                     loggingFunctions,
                     genomeFunctions, predictionFunctions, subprocess_util,
                     __version__)


# setup the application logging
LOG = loggingFunctions.create_logger()

def create_temporary_directory(output_directory):
    characters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    random_string = "".join(random.choice(characters) for _ in range(9))
    temp_dir = os.path.join(output_directory,"tmp_"+random_string)
    return temp_dir

def check_database_struct(db, dbpath):
    levelOneKeys = ["version","date","O","H"]

    for key in levelOneKeys:
        if key not in db.keys():
            raise ValueError("Required key \'{}\' not present in database. Check database at {}".format(key,dbpath))
    if len(db["O"].keys()) == 0:
        raise ValueError("No O antigen alleles found in database. Check database at {}".format(dbpath))
    if len(db["O"].keys()) == 0:
        raise ValueError("No O antigen alleles found in database. Check database at {}".format(dbpath))

    levelTwoKeys = ["gene","desc","allele","seq","MinPident","MinPcov"]

    for antigen in ["O","H"]:
        allele = list(db[antigen].keys())[0]
        for key in levelTwoKeys:
            if key not in db[antigen][allele].keys():
                raise ValueError("Required key \'{}\' not present in database for antigen {} allele {}."
                                 "Check database at {}.".format(key, antigen,allele, dbpath))
    LOG.info("Database structure QC is OK at {}".format(dbpath))

def decompress_gunzip_files(raw_genome_files, temp_dir):
    for idx, g in enumerate(raw_genome_files):
        if 'gz' in g:
            LOG.info(f"Decompression of the gunzip {g} file started ...")
            with open(g, 'rb') as inf, open(f'{temp_dir}/{os.path.basename(g)[:-3]}', 'w', encoding='utf8') as tof:
                decom_str = gzip.decompress(inf.read()).decode('utf-8')
                tof.write(decom_str)    
            LOG.info(f"Wrote decompressed file to {temp_dir}/{os.path.basename(g)[:-3]}")
            raw_genome_files[idx]=os.path.join(temp_dir,os.path.basename(g)[:-3])
    return raw_genome_files

def run_program():
    """
    Main function for E. coli serotyping.
    Creates all required files and controls function execution.
    :return: success or failure
    """
    LOG.setLevel(logging.INFO)
    args = commandLineOptions.parse_command_line()
    
    
    output_directory = create_output_directory(args)
    
    # Create a file handler for log messages in the output directory for the root thread
    fh = logging.FileHandler(os.path.join(output_directory, 'ectyper.log'), 'w', 'utf-8')
    
    if args.debug:
        fh.setLevel(logging.DEBUG)
        LOG.setLevel(logging.DEBUG)
    else:
        fh.setLevel(logging.INFO)
        
    LOG.addHandler(fh)

    #try to load database
    if args.dbpath:
        dbpath = args.dbpath
    else:
        dbpath = definitions.SEROTYPE_ALLELE_JSON

    if os.path.exists(dbpath) == False:
        LOG.critical("ECTyper Database path {} not found".format(args.dbpath))
        raise ValueError("Path to allele database not found. Path specified {}".format(dbpath))

    try:
        with open(file=dbpath) as fp:
            ectyperdb_dict = json.load(fp)
    except:
        raise Exception("Could not load database JSON file")


    check_database_struct(ectyperdb_dict, dbpath)
    
    with open(definitions.PATHOTYPE_ALLELE_JSON) as fp:
        pathotype_db = json.load(fp)

    LOG.info("Starting ectyper v{} running on O and H antigen allele database v{} ({}) and pathotype database v{}".format(
        __version__, 
        ectyperdb_dict["version"], 
        ectyperdb_dict["date"],
        pathotype_db['version']))
    LOG.debug("Command-line arguments were:\n{}".format(args))
    LOG.info("Output_directory is {}".format(output_directory))
    LOG.info("Command-line arguments {}".format(args))

    # Init MASH species database for species identification
    if speciesIdentification.get_species_mash(args.reference) == False:
        LOG.critical("MASH RefSeq sketch does not exists and was not able to be downloaded. Aborting run ...")
        exit("No MASH RefSeq sketch file found at {}".format(args.reference))    

    # Initialize ectyper temporary directory. If --debug is specified then temp folder will be not be deleted. 
    # Python 3.12 introduced delete = False/True option in tempfile lib, so using explicit code supporting Python < 3.12
    temp_dir = create_temporary_directory(output_directory)
    os.makedirs(temp_dir, exist_ok=True)
   
    LOG.info("Gathering genome files list ...")
    
    input_files_list = genomeFunctions.get_files_as_list(args.input, args.maxdirdepth)
    raw_genome_files = decompress_gunzip_files(input_files_list, temp_dir)
  
    LOG.info(f"Identifying genome file types on {len(raw_genome_files)} inputs ...")
    raw_files_dict = genomeFunctions.identify_raw_files(raw_genome_files,
                                                            args)


    alleles_fasta_file = create_alleles_fasta_file(temp_dir, ectyperdb_dict) #from JSON ectyper database file

    combined_fasta = \
            genomeFunctions.create_combined_alleles_and_markers_file(
                alleles_fasta_file, temp_dir, args.pathotype) #create a fasta reference from O-H alleles and optionally from pathotypes alleles database

    
    bowtie_base = genomeFunctions.create_bowtie_base(temp_dir,
                                                         combined_fasta, args.cores) if \
                                                         raw_files_dict['fastq'] else None #only run this function on raw read inputs
        

    # Assemble any fastq files, get final fasta list
    LOG.info("Assembling final list of fasta files")
    all_fastafastq_files_dict = genomeFunctions.assemble_fastq(raw_files_dict,
                                                         temp_dir,
                                                         combined_fasta,
                                                         bowtie_base,
                                                         args)

    # Verify we have at least one fasta file. Optionally species ID.
    # Get a tuple of ecoli and other genomes
    (ecoli_genomes_dict, other_genomes_dict, filesnotfound_dict) = speciesIdentification.verify_ecoli_and_inputs(
            all_fastafastq_files_dict,
            raw_files_dict['other'],
            raw_files_dict['filesnotfound'],
            args)
    
    
    LOG.info("Standardizing the E.coli genome headers based on file names")
    
    predictions_dict={}; predictions_pathotype_dict={}
    if ecoli_genomes_dict:
        ecoli_genomes_dict = genomeFunctions.get_genome_names_from_files(
                ecoli_genomes_dict,
                temp_dir,
                args
                )
        # Run main serotype prediction function
        predictions_dict = run_prediction(ecoli_genomes_dict,
                                          args,
                                          alleles_fasta_file,
                                          temp_dir,
                                          ectyperdb_dict)
           
          
    if other_genomes_dict:
        other_genomes_dict = genomeFunctions.get_genome_names_from_files(
                other_genomes_dict,
                temp_dir,
                args
                )  
            
                

    # Run pathotype predictions if requested irrespective of --verify option
    if args.pathotype:
        predictions_pathotype_dict = predictionFunctions.predict_pathotype_and_shiga_toxin_subtype(ecoli_genomes_dict, other_genomes_dict,
                                                                                   temp_dir,
                                                                                   args.verify,
                                                                                   args.percentIdentityPathotype, 
                                                                                   args.percentCoveragePathotype, 
                                                                                   args.output, args.debug, pathotype_db)

    # Add empty rows for genomes without a blast result or non-E.coli samples that did not undergo typing
    final_predictions = predictionFunctions.add_non_predicted(
        raw_genome_files, predictions_dict, other_genomes_dict, filesnotfound_dict, ecoli_genomes_dict)
    
    
    for sample in final_predictions.keys():
        final_predictions[sample]["database"] = "v"+ectyperdb_dict["version"] + " (" + ectyperdb_dict["date"] + ")"
        if args.pathotype:
            if sample in predictions_pathotype_dict: #could happen that file is non-fasta and pathotyping is called
                final_predictions[sample]["pathotype"] = "/".join(sorted(predictions_pathotype_dict[sample]['pathotype']))
                for field_name in  [f for f in definitions.PATHOTYPE_TOXIN_FIELDS if 'pathotype_' in f]:
                    final_predictions[sample][field_name] = predictions_pathotype_dict[sample][field_name]
                
                for field_name in [f for f in definitions.PATHOTYPE_TOXIN_FIELDS if 'stx_' in f]:
                    final_predictions[sample][field_name] = predictions_pathotype_dict[sample][field_name]
                    
        if 'O' in final_predictions[sample]: #not all samples will have O-antigen dictionary
            highsimilar_Ogroup = getOantigenHighSimilarGroup(final_predictions,sample)


            if highsimilar_Ogroup:
                final_predictions[sample]['O']['highlysimilargroup'] = highsimilar_Ogroup
                final_predictions[sample]['O']['highlysimilarantigens'] = definitions.OSEROTYPE_GROUPS_DICT[highsimilar_Ogroup]
                final_predictions[sample]['error'] = final_predictions[sample]['error'] + \
                                                         "High similarity O-antigen group {}:{} as per A.Iguchi et.al (PMID: 25428893)".format(
                                                          final_predictions[sample]['O']['highlysimilargroup'],
                                                          "/".join(final_predictions[sample]['O']['highlysimilarantigens']))
            else:
                final_predictions[sample]['O']['highlysimilargroup'] = ''
                final_predictions[sample]['O']['highlysimilarantigens'] = []

        if args.verify:
            final_predictions[sample]["QC"] = predictionFunctions.getQuality_control_results(sample,final_predictions,ectyperdb_dict)
  
    # Store most recent result in working directory
    LOG.info("Reporting final results to output.tsv file ...")
    predictionFunctions.report_result(final_predictions, output_directory,
                                          os.path.join(output_directory,
                                                       'output.tsv'),args)
    
    if args.debug == False:
        shutil.rmtree(temp_dir, ignore_errors=True)
    LOG.info(f"ECTyper has finished successfully. Results available at {os.path.abspath(args.output)}")

def getOantigenHighSimilarGroup(final_predictions, sample):
    pred_Otypes = final_predictions[sample]['O']["serogroup"].split("/") #if call is a mixed call

    highlysimilar_Ogroups = set()
    for group_number in definitions.OSEROTYPE_GROUPS_DICT.keys():
        for Otype in pred_Otypes:
            if Otype in definitions.OSEROTYPE_GROUPS_DICT[group_number]:
                highlysimilar_Ogroups.add(group_number)

    if len(highlysimilar_Ogroups) == 1:
        return next(iter(highlysimilar_Ogroups))
    elif len(highlysimilar_Ogroups) >= 2:
        LOG.warning("O-types belong to different high similarity O-antigen groups: {}".format(highlysimilar_Ogroups))
        return  None
    else:
        return None



def create_output_directory(args):
    """
    Create the output directory for ectyper if does not exist already

    :param output_dir: The user-specified output directory, if any
    :return: The output directory
    """
    # If no output directory is specified for the run, create a one based on time
    


    if args.output is None:
        date_dir = ''.join([
            'ectyper_',
            str(datetime.datetime.now().date()),
            '_',
            str(datetime.datetime.now().time()).replace(':', '.')
        ])
        out_dir = os.path.join(definitions.WORKPLACE_DIR, date_dir)
        LOG.info(f"No output folder specified .... All output will be saved in {out_dir}")
        args.output = out_dir
    else:
        if os.path.isabs(args.output):
            out_dir = args.output 
        else:
            out_dir = os.path.join(definitions.WORKPLACE_DIR, args.output)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # clean previous ECTyper output files if the directory was used in previous runs 
    for file in definitions.OUTPUT_FILES_LIST:
        path2file = os.path.join(out_dir,file)
        if os.path.exists(path2file):
            LOG.info(f"Cleaning ECTyper previous files. Removing previously generated {path2file} ...")
            os.remove(path2file) 
    return out_dir



def create_alleles_fasta_file(temp_dir, ectyperdb_dict):
    """
    Every run, re-create the fasta file of alleles to ensure a single
    source of truth for the ectyper data -- the JSON file.

    :temp_dir: temporary directory for length of program run
    :dbpath: path to O and H antigen alleles database in JSON format
    :return: the filepath for alleles.fasta
    """
    output_file = os.path.join(temp_dir, 'alleles.fasta')

    with open(output_file, 'w') as ofh:
        for a in ["O", "H"]:
            for k in ectyperdb_dict[a].keys():
                ofh.write(">" + k + "\n")
                ofh.write(ectyperdb_dict[a][k]["seq"] + "\n")

    LOG.debug(output_file)
    return output_file


def run_prediction(genome_files_dict, args, alleles_fasta, temp_dir, ectyperdb_dict):
    """
    Serotype prediction (O:H antigens) of all the input files, which have now been properly
    converted to fasta if required, and their headers standardized

    :param genome_files: List of genome files in fasta format
    :param args: program arguments from the commandline
    :param alleles_fasta: fasta format file of the ectyper O- and H-alleles
    :param predictions_file: the output file to store the predictions
    :param temp_dir: ectyper run temp dir
    :return: predictions_dict
    """

    genome_files = [genome_files_dict[k]["modheaderfile"] for k in genome_files_dict.keys()]
    # Divide genome files into groups and create databases for each set
    per_core = int(len(genome_files) / args.cores) + 1
    group_size = 50 if per_core > 50 else per_core

    genome_groups = [
        genome_files[i:i + group_size]
        for i in range(0, len(genome_files), group_size)
    ]
    LOG.debug(f"Using O- and H- antigen allele database FASTA file for serotype predictions located at {alleles_fasta}")
    gp = partial(genome_group_prediction, alleles_fasta=alleles_fasta,
                 args=args, temp_dir=temp_dir, ectyperdb_dict=ectyperdb_dict)

    predictions_dict = {}
    with Pool(processes=args.cores) as pool:
        results = pool.map(gp, genome_groups)

        # merge the database predictions with the final predictions dict
        for r in results:
            predictions_dict = {**r, **predictions_dict}
           
    for genome_name in predictions_dict.keys():
        predictions_dict[genome_name]["species"] = "-"
        predictions_dict[genome_name]["species_mash_hash_ratio2ref"] = "-"
        predictions_dict[genome_name]["species_mash_dist2ref"] = "-"
        predictions_dict[genome_name]["species_mash_top_reference"] = "-"
        try:
            predictions_dict[genome_name]["species"] = genome_files_dict[genome_name]["species"]
            predictions_dict[genome_name]["species_mash_hash_ratio2ref"] = genome_files_dict[genome_name]["species_mash_hash_ratio2ref"]
            predictions_dict[genome_name]["species_mash_dist2ref"] = genome_files_dict[genome_name]["species_mash_dist2ref"]
            predictions_dict[genome_name]["species_mash_top_reference"] = genome_files_dict[genome_name]["species_mash_top_reference"]
            predictions_dict[genome_name]["error"] = genome_files_dict[genome_name]["error"]
        except KeyError as e:
            predictions_dict[genome_name]["error"] = "Error: "+str(e)+" in "+genome_name
            LOG.error(f"Failed on {genome_name} sample that does not exist in the 'genome_files_dict' dictionary with the {e} error")

    return predictions_dict


def genome_group_prediction(g_group, alleles_fasta, args, temp_dir, ectyperdb_dict):
    """
    For each genome group, run BLAST, get results, filter and make serotype predictions
    :param g_group: The group of genomes being analyzed
    :param alleles_fasta: fasta format file of the ectyper O- and H-alleles
    :param args: commandline args
    :param temp_dir: ectyper run temp dir
    :return: dictionary of the results for the g_group
    """
    
    # create a temp dir for blastdb -- each process gets its own directory
    LOG.setLevel(logging.INFO) #set level to info as by default only WARNING level is set at init time
    temp_dir_group = create_temporary_directory(temp_dir)
    LOG.info("Creating blast O- and H- antigen database from {}".format(g_group))
    blast_db = os.path.join(temp_dir_group, "blastdb_")
    blast_db_cmd = [
            "makeblastdb",
            "-in", ' '.join(g_group),
            "-dbtype", "nucl",
            "-title", "ectyper_blastdb",
            "-out", blast_db]
    subprocess_util.run_subprocess(blast_db_cmd)

    LOG.info("Starting blast alignment on O- and H- antigen database {}".format(g_group))
    blast_output_file = blast_db + ".output"
    bcline = [
            'blastn',
            '-query', alleles_fasta,
            '-db', blast_db,
            '-out', blast_output_file,
            '-max_hsps', "1",
            '-outfmt',
            "6 qseqid qlen sseqid length pident sstart send sframe qcovhsp bitscore sseq",
            '-word_size', "11"
        ]
    subprocess_util.run_subprocess(bcline)

    LOG.debug("Starting serotype prediction for database {}".format(g_group))
    LOG.debug("BLAST results output file {}".format(blast_output_file))




    db_prediction_dict, blast_output_df = predictionFunctions.predict_serotype(
            blast_output_file,
            ectyperdb_dict,
            args);
    
    blast_output_file_path = os.path.join(args.output,f"blastn_output_alleles.txt")
    if os.path.exists(blast_output_file_path) == False:
        blast_output_df[sorted(blast_output_df.columns)].to_csv(blast_output_file_path , sep="\t", index=False)
        LOG.info("BLAST output file against reference alleles is written at {}".format(blast_output_file_path))
    else:
        blast_output_df[sorted(blast_output_df.columns)].to_csv(blast_output_file_path , mode="a", header=False, sep="\t", index=False)
        LOG.info("Appending BLAST output file against reference alleles at {}".format(blast_output_file_path))
    

    return db_prediction_dict
