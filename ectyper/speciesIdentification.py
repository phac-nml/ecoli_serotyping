import logging
import os
import tempfile
from ectyper import genomeFunctions, definitions, subprocess_util
import re

LOG = logging.getLogger(__name__)


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

def is_shigella(f, args):
    """
    Checks whether the given genome is Shigella or not
    :param f: path to the genome file (FASTA,FASTQ)
    :param args: command line arguments passed to ecTyper during execution
    :return: Boolean True or False
    """
    speciesname = get_species(f,args)
    if re.match("shigella",speciesname,re.IGNORECASE):
        return True
    else:
        return False

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
        subprocess_util.run_subprocess(bcline)

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

    head_cmd = [
        'head',
        '-n', '1'
    ]

    head_output = subprocess_util.run_subprocess(head_cmd,
                                                 input_data=sort_output.stdout)
    top_match = head_output.stdout.decode("utf-8").split()[0]
    LOG.info(top_match)

    m = re.match(r"(GCF_\d+)", top_match)
    refseq_key = None
    if m:
        refseq_key = m.group(1)
    else:
        LOG.critical("Unknown key from MASH species search")
        exit("MASH error")

    LOG.info(refseq_key)
    grep_cmd = [
        'grep',
        refseq_key,
        definitions.REFSEQ_SUMMARY
    ]
    grep_output = subprocess_util.run_subprocess(grep_cmd)

    species = grep_output.stdout.decode("utf-8").split('\t')[7]
    LOG.info(species)

    return species

def getSampleName(file):
    # get only the name of the file for use in the fasta header
    file_base_name = os.path.basename(file)
    file_path_name = os.path.splitext(file_base_name)[0]
    n_name = file_path_name.replace(' ', '_')  # sample name
    return n_name

def verify_ecoli(fasta_files, ofiles, args, temp_dir):
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

    for f in fasta_files:
        sampleName = getSampleName(f)
        if args.verify:
            if is_ecoli(f, temp_dir):
                #ecoli_files.append(f)
                ecoli_files_dict[sampleName] = {"species":"Escherichia coli","filepath":f}
            elif is_shigella(f, temp_dir):
                other_files_dict[sampleName] = {"species": "Shigella","filepath":f}
            else:
                if args.refseq:
                     speciesname = get_species(f, args)
                     other_files_dict[sampleName] = {"species":speciesname,"filepath":f}
                else:
                    other_files_dict[sampleName] = {"species":"Failed E. coli species confirmation based on 10 E.coli specific markers"}
        else:
            #ecoli_files.append(f)
            ecoli_files_dict[f] = {"species":get_species(f, args),"filepath":f}

    for bf in ofiles:
        sampleName = getSampleName(bf)
        other_files_dict[sampleName] = {"message":"Non fasta / fastq file","filepath":bf}

    return ecoli_files_dict, other_files_dict
