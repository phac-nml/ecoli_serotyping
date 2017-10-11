import logging
import multiprocessing
import os

from ectyper import blastFunctions, definitions, subprocess_util

log = logging.getLogger(__name__)

def is_ecoli_genome(fasta_file, args, iden_file=None):
    '''
    Return True if file is classified as ecoli by ecoli markers, otherwise False

    Args:
        fasta_file (str): path to valid fasta genome file
        args (arguments): console arguments
        iden_file (str): Optional path to valid fasta/fastq file
            Used for species identification

    Returns:
        bool: output
    '''
    if iden_file is None:
        iden_file = fasta_file
    num_hit = get_num_hits(fasta_file, args)
    if num_hit < 3:
        log.info("%s is not a valid e.coli genome file", fasta_file)
        if args.species:
            species = get_species(iden_file)
            log.info("%s is identified as genome of %s", fasta_file, species)
        return False
    log.info("%s is a valid e.coli genome file", fasta_file)
    return True

def get_num_hits(target, args):
    '''
    Return number of matching hits when query the reference genome
        on the target genome
    
    Args:
        target(str): target genome
        args (arguments): console arguments

    Returns:
        int: number of hits found
    '''
    num_hit = 0
    try:
        blast_db = blastFunctions.create_blast_db([target])
        result = blastFunctions.run_blast_for_identification(
            definitions.ECOLI_MARKERS,
            blast_db
        )
        with open(result) as handler:
            log.debug("get_num_hits() output:")
            for line in handler:
                log.debug(line)
                num_hit += 1
        log.info("%s aligned to %d marker sequences", target, num_hit)
    except SystemExit:
        pass
    return num_hit

def get_species(file):
    '''
    Given a fasta/fastq file, return the most likely species identification
    
    Args:
        file(str): fasta/fastq file input

    Returns:
        str: name of estimated species
    '''
    if not os.path.isfile(definitions.REFSEQ_SKETCH):
        log.info("No refseq found." +
            "Download refseq from " +
            "https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh" +
            " then put it in ectyper/Data/"
        )
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
    except IndexError:
        log.warning('No matching species found with distance estimation.')
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
            print(screen_output.split('\t'))
            species = screen_output.split('\t')[5]
        except IndexError:
            log.warning('No matching species found with distance screening either.')
    return species
