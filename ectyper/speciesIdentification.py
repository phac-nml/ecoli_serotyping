import logging
import os
import tempfile
from ectyper import definitions, subprocess_util
import re
import requests
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
from Bio import SeqIO
import time #for file age calculations


LOG = logging.getLogger(__name__)


def bool_downloadMashSketch(targetpath):
    if os.path.exists(targetpath) == False or os.path.exists(targetpath+'.txt') == False :
       return True
    else:
        return False

def setLockFile(lockfilepath):
    if os.path.exists(lockfilepath) == False:
        try:
            open(file=lockfilepath, mode="w").close()
            LOG.info("Placed lock file at {}".format(lockfilepath))
        except Exception as e:
            LOG.error("Failed to place a lock file at {}. Database diretory can not be accessed. Wrong path?".format(
                lockfilepath))
            LOG.error("{}".format(e))
            raise FileNotFoundError("Failed to place a lock file at {}".format(lockfilepath))
    else:
        while os.path.exists(lockfilepath):
            elapsed_time = time.time() - os.path.getmtime(lockfilepath) #get file age
            LOG.info("Lock file found at {}. Waiting for other processes to finish database init ...".format(lockfilepath))
            LOG.info("Elapsed time {} min. Will continue processing after 10 min mark even if lock file is present.".format(int(elapsed_time / 60)))
            if elapsed_time >= 600:
                LOG.info("Elapsed time {} min. Assuming previous process completed all init steps. Continue ...".format(int(elapsed_time / 60)))
                if os.path.exists(lockfilepath):  # if previous process failed and 10 min passed since, remove the lock
                    os.remove(lockfilepath)
                break
            time.sleep(60)  # recheck every 1 min if lock file was removed
        LOG.info("Lock file doest not exist or is removed by other process. Continue with databases download ...")

def get_species_mash(targetpath):
    """
    Get MASH sketch of genomes for species identification and check that the most recent version is installed
    :return returns boolean value depending on success of failure to download the MASH sketch
    """


    lockfilepath = os.path.join(os.path.dirname(__file__),"Data/.lock") #prevents race condition
    if os.path.exists(lockfilepath):
        os.remove(lockfilepath)


    if bool_downloadMashSketch(targetpath):
        LOG.info("MASH species id sketch is missing and needs to be downloaded ...")
        setLockFile(lockfilepath)
        for url in definitions.MASH_URLS:
            try:
                if os.path.exists(targetpath) == False:
                    LOG.info("Downloading ~900MB from {}.".format(url))
                    response = requests.get(url,timeout=10, verify=False)
                    response.raise_for_status()
                    LOG.info("Response code for MASH sketch download for mirror {} is {}".format(url, response.status_code))
                    if response.status_code == 200:
                        with open(file=targetpath, mode="wb") as fp:
                            fp.write(response.content)
                generate_sketch_info_summary(targetpath)
                os.remove(lockfilepath)
            except Exception as e:
                LOG.error("Failed to download {}.\nError msg {}".format(url,e))
                pass

            #checks if download was successful and of the right size
            if bool_downloadMashSketch(targetpath) == False:
                LOG.info("Successfully downloaded MASH sketch from {} for species verification".format(url))
                return True
            else:
                LOG.error("Something went wrong with the file download from {}. Trying next mirror ...".format(url))
        LOG.error("All mirrors tried but none worked. Check you Internet connection or download manually (see README)...")
        os.remove(lockfilepath) # remove lock file
        return False #if all mirrors failed

    else:
        LOG.info("MASH species id sketch is in good health and does not need to be downloaded".format(
            targetpath
        ))
        return True

def generate_sketch_info_summary(mash_sketch_path):
    mash_sketch_metadata_file = mash_sketch_path+'.txt'
    process = subprocess_util.run_subprocess(['mash','info', '-t', mash_sketch_path])
    if process.returncode == 0:
        with open(mash_sketch_path+".txt", "w") as fp:
            fp.write(process.stdout.decode('utf-8'))
        LOG.info(f'Successfully generated sketch metadata file at {mash_sketch_metadata_file}')



def is_escherichia_genus(speciesname):
    if re.match(r"Escherichia",speciesname):
        return  True
    else:
        return  False

def get_num_hits(target, temp_dir):
    """
    Identify the number of E. coli specific markers carried by the target
    genome.
    :param target: The genome file under analysis
    :param temp_dir: ectyper run temp_dir
    :return: The number of E. coli specific markers found
    """

    num_hit = 0
    name = os.path.splitext(os.path.split(target)[1])[0]

    with tempfile.TemporaryDirectory(dir=temp_dir) as tdir:
        blast_db = os.path.join(tdir, name)
        blast_db_cmd = [
            "makeblastdb",
            "-in", target,
            "-dbtype", "nucl",
            "-title", "ecoli_test",
            "-out", blast_db]
        subprocess_util.run_subprocess(blast_db_cmd)

        result_file = blast_db + ".output"
        bcline = [
            'blastn',
            '-query', definitions.ECOLI_MARKERS,
            '-db', blast_db,
            '-out', result_file,
            '-perc_identity', "90",
            '-qcov_hsp_perc', "90",
            '-max_target_seqs', "1",
            '-outfmt', "6 qseqid qlen sseqid length pident sstart send sframe",
            '-word_size', "11"
        ]
        subprocess_util.run_subprocess(bcline) #run blast against reference database

        with open(result_file) as handler:
            LOG.debug("get_num_hits() output:")
            LOG.debug("qseqid\tqlen\tsseqid\tlength\tpident\tsstart\tsend\tsframe")
            for line in handler:
                LOG.debug(line.strip())
                num_hit += 1
        LOG.debug("{0} aligned to {1} marker sequences".format(target, num_hit))

    return num_hit


def get_species(file, args, cores=1):
    """
    Given a fasta/fastq file, return the most likely species identification

    Args:
        file (str): fasta/fastq file input

    Returns:
        str: name of estimated species
    """

    top_match="-"; top_match_dist="-"; top_match_hashratio="-"; species="-"
    sketch_metadata_file = args.reference+'.txt' 
    if os.path.exists(sketch_metadata_file) == False:
        raise FileNotFoundError(f'Missing required species ID sketch at {sketch_metadata_file} path')

    mash_cmd = [
        'mash', 'dist', '-p', f"{cores}",
        args.reference,
        file
    ]
    mash_output = subprocess_util.run_subprocess(mash_cmd)

    sort_cmd = [
        'sort',
        '-gk3'
    ]
    sort_output = subprocess_util.run_subprocess(sort_cmd,
                                                 input_data=mash_output.stdout)
    

    if args.debug:
        LOG.debug("Wrote MASH against reference sketch results to {}".format(args.output))
        with open(file=args.output+"/mash_output.txt", mode="w") as fp:
            fp.write(sort_output.stdout.decode("utf-8"))
        fp.close()


    head_cmd = [
        'head',
        '-n', '5'
    ]

    head_output = subprocess_util.run_subprocess(head_cmd,
                                                 input_data=sort_output.stdout)
    top_hit_lines = head_output.stdout.decode("utf-8").split('\n')
   


    if len(top_hit_lines) < 1:
        LOG.warning('For {file} no hits returned by MASH species id sketch search. Species identification failed!')
    else:
        LOG.info('For {} following top hits and hash ratios returned by MASH {}'.format(file,
                                                                        [(top_hit_line.split("\t")[0],top_hit_line.split("\t")[4]) for top_hit_line in top_hit_lines if len(top_hit_line.split("\t")[0])>0]))


    for top_hit_line in top_hit_lines:

        top_hit_line_elements = top_hit_line.split()

        if len(top_hit_line_elements) < 5:
            LOG.warning("No columns in the mash results output to split. Species identification failed!")
            continue

        top_match = top_hit_line_elements[0]; top_match_dist = top_hit_line_elements[2]; top_match_hashratio = top_hit_line_elements[4]
        matched_hashes = top_match_hashratio.split('/')[0]
        matched_meta_line = subprocess_util.run_subprocess(['grep',top_match, sketch_metadata_file],
                                                      ignorereturncode=True).stdout.decode('utf-8').split('\t')
        if len(matched_meta_line) == 4 and matched_hashes != '0':
            m=re.search('s__(.+)',matched_meta_line[3])
            if m:
               species = m.group(1).strip('"')
               LOG.info(
        "MASH species top hit {} identified as {} with distance {} to {} and shared hashes ratio {}".format(top_match, species, top_match_dist, file,
                                                                                            top_match_hashratio))
               LOG.info("MASH dist predicted species name: '{}' based on species ID sketch {}".format(species, args.reference))
        else:
            LOG.warning(f"Could not determine species based on MASH distance for {file}")
            species = "-"
        return species    
        

def getSampleName(file):
    # get only the name of the file for use in the fasta header
    file_base_name = os.path.basename(file)
    file_path_name = os.path.splitext(file_base_name)[0]
    n_name = file_path_name.replace(' ', '_')  # sample name
    return n_name

def is_valid_fasta_file(fasta, sampleName):
    # try to read the first sequence of FASTA file and make a format validity decision. No reason to check all reads
    for contig in SeqIO.parse(fasta, "fasta").records:
        if contig.seq != '':
            LOG.debug(f'{sampleName}: input file {fasta} is a valid FASTA')       
            return True
        else:
            LOG.warning(f'{sampleName}: input FASTA file {fasta} format is invalid FASTA. Skipping further analyses ...')
            return False          
            
    
def verify_ecoli_and_inputs(fasta_fastq_files_dict, ofiles, filesnotfound, args):
    """
    Verifying the E. coli-ness of the genome files and validity of file inputs
    :param fasta_files: [] of all fasta files
    :param ofiles: [] of all non-fasta files
    :param args: Command line arguments
    :return: ecoli_genomes dictionary, other_genomes dictionary, notfound_files dictionary
    """
    LOG.info("Verifying the E.coli-ness of the genome files and validity of file inputs")
    ecoli_files_dict = {}
    other_files_dict = {}
    filesnotfound_dict = {}

    fasta_files = fasta_fastq_files_dict.keys()
    for fasta in fasta_files:
        sampleName = getSampleName(fasta)
        speciesname = "-"

        if is_valid_fasta_file(fasta, sampleName) == False:
            failverifyerrormessage = f"Sample {sampleName} FASTA file ({fasta}) is empty. This could happen when FASTA file generated from FASTQ input lacks raw reads mapping to O- and H- antigens database or input FASTA is empty/corrupted. Please check sequence input file of {sampleName}"
            other_files_dict[sampleName] = {"species":speciesname,"filepath":fasta,"error":failverifyerrormessage}
            return ecoli_files_dict, other_files_dict, filesnotfound_dict

        if sampleName in ecoli_files_dict or sampleName in other_files_dict:
            error_msg = "Duplicated parsed filenames found ('{}'). Offending file paths {}. Only unique file names are supported in batch mode".format(
                                sampleName, [file for file in fasta_files if sampleName in file]
            )
            LOG.error(error_msg)
            raise ValueError(error_msg)

        #do species always regardless of --verify param. Do prediction on fastq files if available for better accuracy
        if fasta_fastq_files_dict[fasta]:
            fastq_file = fasta_fastq_files_dict[fasta]
            speciesname = get_species(fastq_file, args, args.cores)
        else:
            speciesname = get_species(fasta, args, args.cores)

        if args.verify:
            failverifyerrormessage = "Sample identified as " + speciesname + ": serotyping results are only available for E.coli samples." \
                                                                             "If sure that sample is E.coli run without --verify parameter."
            if re.match("Escherichia coli", speciesname):
                ecoli_files_dict[sampleName] = {"species":speciesname,
                                                "filepath":fasta, "error": ""}
            elif is_escherichia_genus(speciesname):
                other_files_dict[sampleName] = {"species":speciesname,"filepath":fasta,"error":failverifyerrormessage}
            else:
                other_files_dict[sampleName] = {"species":speciesname, "filepath":fasta, "error":failverifyerrormessage}
        else:
            ecoli_files_dict[sampleName] = {"species": speciesname,
                                            "filepath": fasta, "error": ""}

    for bf in ofiles:
        sampleName = getSampleName(bf)
        LOG.warning("Non fasta / fastq file.")
        other_files_dict[sampleName] = {"error":"Non fasta / fastq file. ","filepath":bf,"species":"-"}

    for file in filesnotfound:
        sampleName = getSampleName(file)
        filesnotfound_dict[sampleName]={"error":"File {} not found!".format(file)}

    return ecoli_files_dict, other_files_dict,filesnotfound_dict