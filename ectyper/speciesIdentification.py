import logging
import multiprocessing
import os
import tempfile
from Bio import SeqIO

from ectyper import genomeFunctions, blastFunctions, definitions, subprocess_util

LOG = logging.getLogger(__name__)

def is_ecoli_genome(iden_file, genome_file=None, mash=False):
    '''
    Return True if file is classified as ecoli by ecoli markers, otherwise False

    Args:
        iden_file (str): path to valid fasta genome file
        genome_file (str): Optional path to valid fastq file for reads
        mash (bool): Optional input to decide whether to use mash if genome is
                     identified as non-ecoli

    Returns:
        bool: True if iden_file is ecoli, False otherwise
    '''
    if genome_file is None:
        genome_file = iden_file
    num_hit = get_num_hits(iden_file)
    if num_hit < 3:
        LOG.info(
            "{0} is identified as "
            "an invalid e.coli genome file "
            "by marker approach".format(os.path.basename(iden_file)))
        if mash:
            species = get_species(genome_file)
            LOG.info(
                "{0} is identified as genome of "
                "{1} by mash approach".format(os.path.basename(iden_file), species))
        return False
    LOG.debug("{0} is a valid e.coli genome file".format(os.path.basename(iden_file)))
    return True

def get_num_hits(target):
    '''
    Return number of matching hits when query the reference genome
        on the target genome

    Args:
        target (str): target genome

    Returns:
        int: number of hits found
    '''
    num_hit = 0
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            blast_db = blastFunctions.create_blast_db([target], temp_dir)
            result = blastFunctions.run_blast_for_identification(
                definitions.ECOLI_MARKERS,
                blast_db
            )
            with open(result) as handler:
                LOG.debug("get_num_hits() output:")
                for line in handler:
                    LOG.debug(line)
                    num_hit += 1
            LOG.debug("{0} aligned to {1} marker sequences".format(target, num_hit))
    except SystemExit:
        pass
    return num_hit

def get_species(file):
    '''
    Given a fasta/fastq file, return the most likely species identification

    Args:
        file (str): fasta/fastq file input

    Returns:
        str: name of estimated species
    '''
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
    '''
    Given a fasta/fastq file with one sequence, return the most likely species
    identification

    Args:
        file (str): fasta/fastq file input

    Returns:
        str: name of estimated species
    '''
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
        mash_output = subprocess_util.run_subprocess(' '.join(cmd), is_shell=True)
        ass_acc_num = '_'.join(mash_output.split('\t')[1].split('_')[:2])
        cmd = [
            'grep -E',
            ass_acc_num,
            definitions.REFSEQ_SUMMARY
        ]
        summary_output = subprocess_util.run_subprocess(' '.join(cmd), is_shell=True)
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
            screen_output = subprocess_util.run_subprocess(' '.join(cmd), is_shell=True)
            LOG.debug(screen_output.split('\t'))
            species = screen_output.split('\t')[5].split('\n')[0]
        except Exception as err2:
            LOG.warning(
                'No matching species found with distance screening either:{}'.format(err2)
            )
    return species
