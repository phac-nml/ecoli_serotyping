import logging
import os
import tempfile
from ectyper import genomeFunctions, definitions, subprocess_util
import re
import requests
import time #for file age calculations
import portalocker

LOG = logging.getLogger(__name__)

def get_fileAge(file):
    if os.path.exists(file):
        age = int(time.time() - os.stat(file).st_mtime)
        LOG.info("RefSeq MASH refseq.genomes.k21s1000.msh age is {} days downloaded on {} located {}".format(int(age / 86400),
                                                                                                  time.strftime('%Y-%m-%d at %H:%M:%S', time.localtime(os.path.getmtime(file))), file))
        return age #file age in seconds
    else:
        LOG.error("RefSeq MASH refseq.genomes.k21s1000.msh does not exist. Was not able to get file age.")
        return 0

def bool_downloadMashRefSketch(targetpath):
    if os.path.exists(targetpath) == False:
       return True
    #if the file size is smaller than 700MB, re-download mash sketch from Internet
    elif os.path.getsize(targetpath) < 700000000:
       return True
    #re-download mash sketch after 6 months (16070400 seconds)
    elif get_fileAge(targetpath) > 16070400:#
        LOG.warning("MASH Sketch exists and is >6 month old and needs to be re-downloaded")
        return True
    else:
        return False

def get_refseq_mash():
    """
    Get MASH sketch of refseq genomes for species identification and check that the most recent version is installed
    :return returns boolean value depending on success of failure to download RefSeq MASH sketch
    """

    urls=["https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh",
          "https://share.corefacility.ca/index.php/s/KDhSNQfhE6npIyo/download",
          "https://gitlab.com/kbessonov/ectyper/raw/master/ectyper/Data/refseq.genomes.k21s1000.msh"]

    targetpath = os.path.join(os.path.dirname(__file__),"Data/refseq.genomes.k21s1000.msh")

    if bool_downloadMashRefSketch(targetpath):
        for url in urls:
            LOG.info("Trying to download MASH sketch from {}.".format(url))
            try:
                LOG.info("Downloading ~700MB from {}.".format(url))
                response = requests.get(url,timeout=10, verify=False)
                response.raise_for_status()
                if response.status_code == 200:
                  with portalocker.Lock(filename=targetpath, mode="wb", flags=portalocker.LOCK_EX) as fp:
                      fp.write(response.content); fp.flush()
                  download_assembly_summary()
            except Exception as e:
                LOG.error("Failed to download refseq.genomes.k21s1000.msh from {}.\nError msg {}".format(url,e))
                pass

            # checks if download was successful and of the right size
            if bool_downloadMashRefSketch(targetpath) == False:
                LOG.info("Sucessfully downloaded RefSeq MASH sketch from {} for species verification".format(url))
                return True
            else:
                LOG.error("Something went wrong with the file download or downloaded file is truncated/corrupted from {}. Trying next mirror ...".format(url))
        return False #if all mirrors failed

    else:
        assemblysummarypath = os.path.join(os.path.dirname(__file__), "Data/assembly_summary_refseq.txt")
        if os.path.exists(assemblysummarypath) == False:
            download_assembly_summary()
        LOG.info("RefSeq sketch (refseq.genomes.k21s1000.msh) and assembly meta data (assembly_summary_refseq.txt) is in good health and does not need to be downloaded")
        return True

def download_assembly_summary():
    sourceurls = {"assembly_summary_refseq.txt":"http://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"}

    try:
        for targetfile, url in sourceurls.items():
            targetpath = os.path.join(os.path.dirname(__file__), "Data/"+targetfile)
            response = requests.get(url, timeout=10, verify=False)
            response.raise_for_status()

            if response.status_code == 200:
                with portalocker.Lock(filename=targetpath, mode="w", flags=portalocker.LOCK_EX) as fp:
                    fp.write(response.text); fp.flush()
                LOG.info("Successfully downloaded {} with response code {}".format(targetfile, response.status_code))
            else:
                LOG.critical("Server response error {}. Failed to download {}".format(response.status_code,targetfile))

            if not os.path.exists(targetpath):
                exit(1)
                LOG.critical("The {} file does not exist or is corrupted at {} path".format(targetfile, targetpath))

    except Exception as e:
        print(e)
        LOG.critical("Failed to download or write to a disk assembly_summary_genbank.txt or assembly_summary_refseq.txt.\nError msg {}".format(str(e)))
        pass




def is_ecoli(genome_file, temp_dir,ecolidict):
    """
    Checks whether the given genome is E. coli or not

    :param genome_file: The genome file in fasta format
    :param temp_dir: ectyper run temp_dir
    :return: True or False
    """
    LOG.info("Verifying if sample is a valid E.coli genome based on the 10 E.coli specific markers")
    num_hit = get_num_hits(genome_file, temp_dir)
    ecolidict["nhitsmarkers"] = num_hit

    if num_hit < 9:
        LOG.warning(
            "{0} is identified as an invalid E. coli genome.\nOnly {1} out of 10 E.coli specific markers were identified "
            "by the marker approach of "
            "https://bmcmicrobiol.biomedcentral.com/articles/10.1186/s12866"
            "-016-0680-0#Tab3"
            " where at least 9 out of 10 E. coli specific markers must be "
            "present for a given genome to be considered E.coli".format(os.path.basename(genome_file), num_hit))
        return False
    else:
        LOG.info("Samples is a valid E.coli genome based on the {} markers".format(num_hit))
        return True


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


def get_species(file, args):
    """
    Given a fasta/fastq file, return the most likely species identification

    Args:
        file (str): fasta/fastq file input

    Returns:
        str: name of estimated species
    """

    top_match="-"; top_match_dist="-"; top_match_hashratio="-"; species="-"
    mash_cmd = [
        'mash', 'dist',
        args.refseq,
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
        LOG.debug("Wrote MASH output at {}".format(os.getcwd()))
        with open(file="mash_output.txt", mode="w") as fp:
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
        LOG.warning("No hits returned by mash sketch search. Species identification failed!")
    else:
        LOG.info("Following top hits returned by MASH {}".format([top_hit_line.split("\t")[0] for top_hit_line in top_hit_lines]))


    for top_hit_line in top_hit_lines:

        top_hit_line_elements = top_hit_line.split()

        if len(top_hit_line_elements) < 5:
            LOG.warning("No columns in the mash results ouptut to split. Species identification failed!")
            continue

        top_match = top_hit_line_elements[0]; top_match_dist = top_hit_line_elements[2]; top_match_hashratio = top_hit_line_elements[4]

        m = re.match(r"(GCF_\d+)", top_match)

        top_match_hashratio_tuple = re.findall(r'(\d+)\/(\d+)', top_match_hashratio)[0]
        top_match_sharedhashes = int(top_match_hashratio_tuple[0])

        if m:
            refseq_key = m.group(1)
        else:
            LOG.warning("Could not detemine species based on MASH Distance"
                        "Could not extract GCF_# accession number from the MASH dist results.".format(top_match))
            continue #try other top match

        if m is None or top_match_sharedhashes < 100:
            LOG.warning("\nTop MASH sketch hit {} with {} shared hashes."
                        "\nCould not assign species based on MASH distance to reference sketch file.\n"
                        "Either:\n"
                        "1. MASH sketch meta data accessions do not start with the GCF_ prefix or\n"
                        "2. Number of shared hashes to reference is less than 100 (i.e. too distant).\n"
                        .format(top_match,top_match_hashratio))
            species = "-" # if after top 10 genome IDs still no accession number match, give up
            return species

        LOG.info(refseq_key)
        grep_cmd = [
            'grep',
            refseq_key,
            definitions.REFSEQ_SUMMARY
        ]
        grep_output = subprocess_util.run_subprocess(grep_cmd, ignorereturncode=True)
        grep_output_decoded = grep_output.stdout.decode("utf-8").split('\t')

        if grep_output_decoded and len(grep_output_decoded) >= 8:
            species = grep_output.stdout.decode("utf-8").split('\t')[7]
        else:
            species = "-"

        if not species == "-":
            break  #no need to continue looping if top hit species is found

    LOG.info(
        "MASH species RefSeq top hit {} with distance {} and shared hashes ratio {}".format(top_match, top_match_dist,
                                                                                            top_match_hashratio))
    LOG.info("MASH dist predicted species name: {}".format(species))

    return species

def getSampleName(file):
    # get only the name of the file for use in the fasta header
    file_base_name = os.path.basename(file)
    file_path_name = os.path.splitext(file_base_name)[0]
    n_name = file_path_name.replace(' ', '_')  # sample name
    return n_name

def verify_ecoli(fasta_fastq_files_dict, ofiles, args, temp_dir):
    """
    Verifying the _E. coli_-ness of the genome files
    :param fasta_files: [] of all fasta files
    :param ofiles: [] of all non-fasta files
    :param args: Command line arguments
    :param temp_dir: ectyper run temp_dir
    :return: ([ecoli_genomes], {file:species})
    """

    #ecoli_files = []
    ecoli_files_dict = {}
    other_files_dict = {}

    fasta_files = fasta_fastq_files_dict.keys()
    for fasta in fasta_files:
        sampleName = getSampleName(fasta)
        speciesname = "-"
        numofecolimarkers = {"ecolimarkers": 0} #number of E.coli specific markers

        if args.verify:

            # assign the default RefSeq sketch if none is specified
            if args.refseq == None:
                args.refseq = os.path.join(os.path.dirname(__file__),
                                           "Data/refseq.genomes.k21s1000.msh")

            #do species prediction of fastq files if available for better accuracy
            if fasta_fastq_files_dict[fasta]:
                fastq_file = fasta_fastq_files_dict[fasta]
                speciesname = get_species(fastq_file, args)
            else:
                speciesname = get_species(fasta, args)

            failverifyerrormessage = "Sample identified as " + speciesname + ": serotyping results are only available for E.coli samples"
            if re.match("Escherichia coli", speciesname):
                numofecolimarkers["ecolimarkers"] = get_num_hits(fasta,temp_dir) #get number of E.coli markers
                ecoli_files_dict[sampleName] = {"species":speciesname,"ecolimarkers":numofecolimarkers["ecolimarkers"],
                                                "filepath":fasta, "error": "-"}
            elif is_ecoli(fasta, temp_dir,numofecolimarkers):
                ecoli_files_dict[sampleName] = {"species":"Escherichia coli", "ecolimarkers":numofecolimarkers["ecolimarkers"],
                                                "filepath":fasta,"error":failverifyerrormessage}
            elif is_escherichia_genus(speciesname):
                other_files_dict[sampleName] = {"species":speciesname,"filepath":fasta,"error":failverifyerrormessage}
            else:
                other_files_dict[sampleName] = {"species":speciesname, "error":failverifyerrormessage}
        else:
            ecoli_files_dict[sampleName] = {"species": speciesname,"nmarkers":numofecolimarkers["nmarkers"],
                                            "filepath": fasta, "error": "-"}

    for bf in ofiles:
        sampleName = getSampleName(bf)
        LOG.warning("Non fasta / fastq file")
        other_files_dict[sampleName] = {"error":"Non fasta / fastq file","filepath":bf,"species":"-"}

    return ecoli_files_dict, other_files_dict
