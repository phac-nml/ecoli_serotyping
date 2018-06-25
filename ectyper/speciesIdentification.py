import logging
import multiprocessing
import os
import tempfile
from Bio import SeqIO
from ectyper import genomeFunctions, definitions, subprocess_util

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


def get_species(file):
    """
    Given a fasta/fastq file, return the most likely species identification

    Args:
        file (str): fasta/fastq file input

    Returns:
        str: name of estimated species
    """

    LOG.info("Identifying species for {0}".format(file))
    if not os.path.isfile(definitions.REFSEQ_SKETCH):
        LOG.warning("No refseq found.")
        return None
    species = 'unknown'
    if genomeFunctions.get_valid_format(file) == 'fasta':
        tmp_file = tempfile.NamedTemporaryFile().name
        basename = os.path.basename(file).replace(' ', '_')
        with open(tmp_file, 'w') as new_fh:
            header = '> {0}\n'.format(basename)
            new_fh.write(header)
            for record in SeqIO.parse(file, 'fasta'):
                new_fh.write(str(record.seq))
                new_fh.write('nnnnnnnnnnnnnnnnnnnn')
        try:
            species = get_species_helper(tmp_file)
        except Exception:
            pass
    if genomeFunctions.get_valid_format(file) == 'fastq':
        species = get_species_helper(file)
    return species


def get_species_helper(file):
    """
    Given a fasta/fastq file with one sequence, return the most likely species
    identification

    Args:
        file (str): fasta/fastq file input

    Returns:
        str: name of estimated species
    """

    species = 'unknown'
    cmd = [
        'mash', 'dist',
        file,
        definitions.REFSEQ_SKETCH,
        '|',
        'sort -gk3 -',
        '|',
        'head -1 -'
    ]
    try:
        mash_output = subprocess_util.run_subprocess(' '.join(cmd))
        ass_acc_num = '_'.join(mash_output.split('\t')[1].split('_')[:2])
        cmd = [
            'grep -E',
            ass_acc_num,
            definitions.REFSEQ_SUMMARY
        ]
        summary_output = subprocess_util.run_subprocess(' '.join(cmd))
        species = summary_output.split('\t')[7]
        return species
    except Exception as err:
        LOG.warning('No matching species found with distance estimation:{0}'.format(err))
        try:
            cmd = [
                'mash screen',
                '-w',
                '-p', str(multiprocessing.cpu_count()//2),
                definitions.REFSEQ_SKETCH,
                file,
                '| sort -gr - | head -1 -'
            ]
            screen_output = subprocess_util.run_subprocess(' '.join(cmd))
            LOG.debug(screen_output.split('\t'))
            species = screen_output.split('\t')[5].split('\n')[0]
        except Exception as err2:
            LOG.warning(
                'No matching species found with distance screening either:{}'.format(err2)
            )
    return species


def verify_ecoli(fasta_files, species):
    """
    Verifying the _E. coli_-ness of the genome files
    :param fasta_files: [] of all fasta files
    :return: List of fasta files
    """

    final_files = []
    for f in fasta_files:
        if is_ecoli(f):
            final_files.append(f)
        else:
            if species:
                get_species(f)
                # Report species prediction

    return final_files




