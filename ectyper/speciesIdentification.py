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
        LOG.info("RefSeq MASH refseq.genomes.k21s1000.msh age is {} days downloaded on {}".format(int(age / 86400),
                                                                                                  time.strftime('%Y-%m-%d at %H:%M:%S', time.localtime(os.path.getmtime(file)))))
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
                  portalocker.Lock(targetpath,timeout=3)
                  with open(file=targetpath,mode="wb") as fp:
                      fp.write(response.content)
                  fp.close()

                  download_RefSeq_assembly_summary()
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
            download_RefSeq_assembly_summary()
        LOG.info("RefSeq sketch is in good health and does not need to be downloaded")
        return True

def download_RefSeq_assembly_summary():
    url = "http://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
    targetpath = os.path.join(os.path.dirname(__file__), "Data/assembly_summary_refseq.txt")
    try:
        response = requests.get(url, timeout=10, verify=False)
        response.raise_for_status()


        if response.status_code == 200:
            portalocker.Lock(targetpath, timeout=3)
            with open(file=targetpath,mode="w") as fp:
                fp.write(response.text)
            fp.close()
            LOG.info("Successfully downloaded assembly_summary_refseq.txt".format(response.status_code))
        else:
            LOG.critical("Server response error {}. Failed to download assembly_summary_refseq.txt".format(response.status_code))

    except Exception as e:
        LOG.critical("Failed to download or write to a disk assembly_summary_refseq.txt from {}.\nError msg {}".format(url, str(e)))
        pass

    if (os.path.exists(targetpath)) == False or os.path.getsize(targetpath) < 54000000:
        exit(1)
        LOG.critical("The assembly_summary_refseq.txt file does not exist or is corrupted at {} path".format(targetpath))

def is_ecoli(genome_file, temp_dir):
    """
    Checks whether the given genome is E. coli or not

    :param genome_file: The genome file in fasta format
    :param temp_dir: ectyper run temp_dir
    :return: True or False
    """
    num_hit = get_num_hits(genome_file, temp_dir)

    if num_hit < 3:
        LOG.warning(
            "{0} is identified as an invalid E. coli genome "
            "by the marker approach of "
            "https://bmcmicrobiol.biomedcentral.com/articles/10.1186/s12866"
            "-016-0680-0#Tab3"
            " where at least three E. coli specific markers must be "
            "present".format(os.path.basename(genome_file)))
        return False
    else:
        return True

def is_escherichia_genus(speciesname):
    if re.match("Escherichia",speciesname):
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
            for line in handler:
                LOG.debug(line)
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
        '-n', '1'
    ]

    head_output = subprocess_util.run_subprocess(head_cmd,
                                                 input_data=sort_output.stdout)
    top_hit = head_output.stdout.decode("utf-8").split()
    top_match = top_hit[0]; top_match_dist = top_hit[2]; top_match_hashratio = top_hit[4]

    LOG.info("MASH species RefSeq top hit {} with distance {} and shared hashes ratio {}".format(top_match,top_match_dist,top_match_hashratio))

    m = re.match(r"(GCF_\d+)", top_match)
    #refseq_key = None
    if m:
        refseq_key = m.group(1)
    else:
        LOG.critical("Unknown key from MASH species search")
        return("Undetermined species. Could not map genome accession #{} to species".format(top_match))
        #exit("MASH error")

    LOG.info(refseq_key)
    grep_cmd = [
        'grep',
        refseq_key,
        definitions.REFSEQ_SUMMARY
    ]
    grep_output = subprocess_util.run_subprocess(grep_cmd)

    species = grep_output.stdout.decode("utf-8").split('\t')[7]
    LOG.info("MASH dist predicted species name:{}".format(species))

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

    fasta_files= fasta_fastq_files_dict.keys()
    for fasta in fasta_files:
        sampleName = getSampleName(fasta)

        if args.refseq:
            if fasta_fastq_files_dict[fasta]:
                fastq_file = fasta_fastq_files_dict[fasta]
                speciesname = get_species(fastq_file, args) #would like to submit entire fastq file for more accurate species identification instead of allele-based fasta surrogate
            else:
                speciesname = get_species(fasta, args)
        else:   #if user does not specify a RefSeq sketch, use the default one
            args.refseq = os.path.join(os.path.dirname(__file__),"Data/refseq.genomes.k21s1000.msh") #default sketch
            if fasta_fastq_files_dict[fasta]:
                fastq_file = fasta_fastq_files_dict[fasta]
                speciesname = get_species(fastq_file, args)
            else:
                speciesname = get_species(fasta, args)

        if args.verify:
            if is_ecoli(fasta, temp_dir):
                #ecoli_files.append(f)
                ecoli_files_dict[sampleName] = {"species":"Escherichia coli","filepath":fasta,"error":"-"}
            #elif is_shigella(f, speciesname, args):
            #    other_files_dict[sampleName] = {"species": speciesname, "filepath": f}
            elif is_escherichia_genus(speciesname):
                ecoli_files_dict[sampleName] = {"species":speciesname,"filepath":fasta,"error":"This sample belongs to Escherichia genus so serotyping results might not be accurate"}
            else:
                if args.refseq:
                    other_files_dict[sampleName] = {"species":speciesname,"filepath":fasta,"error":"-"}
                else:
                    other_files_dict[sampleName] = {"species":"Non-Ecoli", "error":"Failed E. coli species confirmation based on 10 E.coli specific markers"}
        else:
            #ecoli_files.append(f)
            ecoli_files_dict[fasta] = {"species":speciesname,"filepath":fasta}

    for bf in ofiles:
        sampleName = getSampleName(bf)
        LOG.warning("Non fasta / fastq file")
        other_files_dict[sampleName] = {"error":"Non fasta / fastq file","filepath":bf,"species":"-"}

    return ecoli_files_dict, other_files_dict
