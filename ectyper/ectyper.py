#!/usr/bin/env python
"""
    Predictive serotyping for _E. coli_.
"""
import os, re
import tempfile
import datetime
import json
import logging
from multiprocessing import Pool
from functools import partial

from ectyper import (commandLineOptions, definitions, speciesIdentification,
                     loggingFunctions,
                     genomeFunctions, predictionFunctions, subprocess_util,
                     __version__)


# setup the application logging

LOG = loggingFunctions.create_logger()

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


def run_program():
    """
    Main function for E. coli serotyping.
    Creates all required files and controls function execution.
    :return: success or failure
    """


    args = commandLineOptions.parse_command_line()
    output_directory = create_output_directory(args.output)

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

    LOG.info("Starting ectyper v{} running on allele database v{} ({})".format(__version__, ectyperdb_dict["version"], ectyperdb_dict["date"]))

    if args.debug:
        LOG.debug("Command-line arguments were:\n{}".format(args))

    # Create a file handler for log messages in the output directory
    fh = logging.FileHandler(os.path.join(output_directory, 'ectyper.log'), 'w', 'utf-8')

    if args.debug:
        fh.setLevel(logging.DEBUG)
        LOG.setLevel(logging.DEBUG)
    else:
        fh.setLevel(logging.INFO)
        LOG.setLevel(logging.INFO)


    LOG.addHandler(fh)
    LOG.info("Output_directory is {}".format(output_directory))
    LOG.info("Command-line arguments {}".format(args))

    # Init RefSeq database for species identification
    if speciesIdentification.get_refseq_mash_and_assembly_summary() == False:
        LOG.critical("MASH RefSeq sketch does not exists and was not able to be downloaded. Aborting run ...")
        exit("No MASH RefSeq sketch file found at {}".format(dbpath))
    else:
        lockfilepath=os.path.join(dbpath,".lock")
        if os.path.exists(lockfilepath):
            os.remove(lockfilepath)

    # Initialize ectyper temporary directory for the scope of this program (it will be deleted when program ends)
    with tempfile.TemporaryDirectory(dir=output_directory) as temp_dir:
        LOG.info("Gathering genome files")
        raw_genome_files = genomeFunctions.get_files_as_list(args.input)

        LOG.info("Identifying genome file types")

        raw_files_dict = genomeFunctions.identify_raw_files(raw_genome_files,
                                                            args)


        alleles_fasta_file = create_alleles_fasta_file(temp_dir, ectyperdb_dict) #from JSON ectyper database file

        combined_fasta = \
            genomeFunctions.create_combined_alleles_and_markers_file(
                alleles_fasta_file, temp_dir)


        bowtie_base = genomeFunctions.create_bowtie_base(temp_dir,
                                                         combined_fasta) if \
                                                         raw_files_dict['fastq'] else None

        # Assemble any fastq files, get final fasta list
        LOG.info("Assembling final list of fasta files")
        all_fastafastq_files_dict = genomeFunctions.assemble_fastq(raw_files_dict,
                                                         temp_dir,
                                                         combined_fasta,
                                                         bowtie_base,
                                                         args)

        # Verify we have at least one fasta file. Optionally species ID.
        # Get a tuple of ecoli and other genomes
        (ecoli_genomes_dict, other_genomes_dict, filesnotfound_dict) = speciesIdentification.verify_ecoli(
            all_fastafastq_files_dict,
            raw_files_dict['other'],
            raw_files_dict['filesnotfound'],
            args,
            temp_dir)


        LOG.info("Standardizing the E.coli genome headers based on file names")
        predictions_dict={}
        if ecoli_genomes_dict:
            ecoli_genomes_dict = genomeFunctions.get_genome_names_from_files(
                ecoli_genomes_dict,
                temp_dir,
                args
                )
            # Run main prediction function
            predictions_dict = run_prediction(ecoli_genomes_dict,
                                          args,
                                          alleles_fasta_file,
                                          temp_dir,
                                          ectyperdb_dict)


        # Add empty rows for genomes without a blast result or non-E.coli samples that did not undergo typing
        final_predictions = predictionFunctions.add_non_predicted(
            raw_genome_files, predictions_dict, other_genomes_dict, filesnotfound_dict, ecoli_genomes_dict)


        for sample in final_predictions.keys():
            final_predictions[sample]["database"] = "v"+ectyperdb_dict["version"] + " (" + ectyperdb_dict["date"] + ")"

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
        LOG.info("Reporting results:")

        predictionFunctions.report_result(final_predictions, output_directory,
                                          os.path.join(output_directory,
                                                       'output.tsv'),args)
        LOG.info("\nECTyper has finished successfully.")

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



def create_output_directory(output_dir):
    """
    Create the output directory for ectyper

    :param output_dir: The user-specified output directory, if any
    :return: The output directory
    """
    # If no output directory is specified for the run, create a one based on
    # time


    if output_dir is None:
        date_dir = ''.join([
            'ectyper_',
            str(datetime.datetime.now().date()),
            '_',
            str(datetime.datetime.now().time()).replace(':', '.')
        ])
        out_dir = os.path.join(definitions.WORKPLACE_DIR, date_dir)
    else:
        if os.path.isabs(output_dir):
            out_dir = output_dir
        else:
            out_dir = os.path.join(definitions.WORKPLACE_DIR, output_dir)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
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
    Serotype prediction of all the input files, which have now been properly
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
    gp = partial(genome_group_prediction, alleles_fasta=alleles_fasta,
                 args=args, temp_dir=temp_dir, ectyperdb_dict=ectyperdb_dict)

    predictions_dict = {}
    with Pool(processes=args.cores) as pool:
        results = pool.map(gp, genome_groups)

        # merge the per-database predictions with the final predictions dict
        for r in results:
            predictions_dict = {**r, **predictions_dict}

    for genome_name in predictions_dict.keys():
        try:
            predictions_dict[genome_name]["species"] = genome_files_dict[genome_name]["species"]
            predictions_dict[genome_name]["error"] = genome_files_dict[genome_name]["error"]
        except KeyError as e:
            predictions_dict[genome_name]["error"] = "Error: "+str(e)+" in "+genome_name

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
    with tempfile.TemporaryDirectory(dir=temp_dir) as temp_dir:
        LOG.debug("Creating blast database from {}".format(g_group))
        blast_db = os.path.join(temp_dir, "blastdb_")
        blast_db_cmd = [
            "makeblastdb",
            "-in", ' '.join(g_group),
            "-dbtype", "nucl",
            "-title", "ectyper_blastdb",
            "-out", blast_db]
        subprocess_util.run_subprocess(blast_db_cmd)

        LOG.debug("Starting blast alignment on database {}".format(g_group))
        blast_output_file = blast_db + ".output"
        bcline = [
            'blastn',
            '-query', alleles_fasta,
            '-db', blast_db,
            '-out', blast_output_file,
           # '-perc_identity', str(args.percentIdentity),
           # '-qcov_hsp_perc', str(args.percentLength),
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
            args)

        blast_output_file_path = args.output+"/blast_output_alleles.txt"
        blast_output_df[sorted(blast_output_df.columns)].to_csv(blast_output_file_path , sep="\t", index=False)
        LOG.info("BLAST output file against reference alleles is written at {}".format(blast_output_file_path))

        return db_prediction_dict
