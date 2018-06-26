import logging
import multiprocessing
import os
import tempfile
from Bio import SeqIO
from ectyper import genomeFunctions, definitions, subprocess_util
import pandas as pd

LOG = logging.getLogger(__name__)


def is_ecoli(genome_file):
    """
    Checks whether the given genome is E. coli or not

    :param genome_file: The genome file in fasta format
    :return: True or False
    """

    num_hit = get_num_hits(genome_file)
    if num_hit < 3:
        LOG.warning(
            "{0} is identified as an invalid E. coli genome"
            "by the marker approach of "
            "https://bmcmicrobiol.biomedcentral.com/articles/10.1186/s12866-016-0680-0#Tab3"
            " where at least three E. coli specific markers must be "
            "present".format(os.path.basename(genome_file)))
        return False
    else:
        return True


def get_num_hits(target):
    """
    Identify the number of E. coli specific markers carried by the target genome.
    :param target: The genome file under analysis
    :return: The number of E. coli specific markers found
    """

    num_hit = 0
    name = os.path.splitext(os.path.split(target)[1])[0]

    with tempfile.TemporaryDirectory() as tdir:
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
    # Download refseq file if species identification is enabled
    refseq_sketch = args.refseq if args.refseq else genomeFunctions.download_refseq()
    species = 'Unknown'

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
    sort_output = subprocess_util.run_subprocess(sort_cmd, input_data=mash_output.stdout)

    head_cmd = [
        'head',
        '-n', '1'
    ]

    head_output = subprocess_util.run_subprocess(head_cmd, input_data=sort_output.stdout)
    top_match = head_output.stdout.decode("utf-8").split()[0]
    LOG.info(top_match)

    return species


def verify_ecoli(fasta_files, args):
    """
    Verifying the _E. coli_-ness of the genome files
    :param fasta_files: [] of all fasta files
    :param args: Command line arguments
    :return: List of fasta files
    """

    final_files = []
    for f in fasta_files:
        if is_ecoli(f):
            final_files.append(f)
        else:
            if args.species:
                get_species(f, args)
                # Report species prediction

    return final_files




