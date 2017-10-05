'''
Genome Utilities
'''
#!/usr/bin/env python

import json
import logging
import os
import re
import tempfile

import Bio
from Bio import SeqIO

from ectyper import definitions
from ectyper import serotypePrediction
from ectyper import subprocess_util
from ectyper import blastFunctions

log = logging.getLogger(__name__)


def get_files_as_list(file_or_directory):
    """
    Creates a list of files from either the given file, or all files within the
    directory specified (where each file name is its absolute path).

    :param file_or_directory: file or directory name given on commandline
    :return files_list: List of all the files found.

    """

    files_list = []
    if file_or_directory == '':
        return files_list

    if os.path.isdir(file_or_directory):
        log.info("Gathering genomes from directory " + file_or_directory)

        # Create a list containing the file names
        for root, dirs, files in os.walk(file_or_directory):
            for filename in files:
                files_list.append(os.path.join(root, filename))

    else:
        log.info("Using genomes in file " + file_or_directory)
        files_list.append(os.path.abspath(file_or_directory))

    return sorted(files_list)


def get_valid_format(file):
    """
    Check using SeqIO if files are valid fasta/fastq format.
    Then return the format.

    Args:
        file (str): path of file
    
    Returns:
        fmt (str): 'fasta', 'fastq', or ''
    """
    file_format = os.path.splitext(file)[1][1:]
    valid_fasta_formats = ['fna', 'fa', 'fasta']
    valid_fastq_formats = ['fq', 'fastq']
    if file_format in valid_fasta_formats:
        file_format = 'fasta'
    elif file_format in valid_fastq_formats:
        file_format = 'fastq'
    else:
        return ''
    for _ in SeqIO.parse(file, file_format):
        log.info("%s is a valid %s format file", file, file_format)
        return file_format
    log.info("%s is not a valid %s format file", file, file_format)
    return ''


def get_genome_names_from_files(files):
    """
    For each file:
    Takes the first header from a fasta file and sends it to the get_genome_name
    function. Returns the name of the genome. If the name of the file is to be
    used as the genome name, creates a temporary file using >lcl|filename as the
    first part of the header.

    :param files: The list of files to get the genome names for
    :return: ([genome names], [file names])
    """

    list_of_genomes = []
    list_of_files = []
    for file in files:
        header = get_fasta_header_from_file(file)
        genome_name = get_genome_name(header)

        # if the header and genome_name are the same, we need to use the filename
        # as the genome name. This means we also need to create a new file adding
        # the filename to each of the headers, so that downstream applications
        # (eg. BLAST) can be used with the filename as genome name.
        if header == '':
            log.fatal('No header for %s', file)
            exit(1)
        if header == genome_name:
            # get only the name of the file for use in the fasta header
            file_path_name = os.path.splitext(os.path.basename(file))
            n_name = file_path_name[0]

            # create a new file for the updated fasta headers
            new_file_tuple = tempfile.mkstemp()
            new_file = new_file_tuple[1]

            # add the new name to the list of files and genomes
            list_of_files.append(new_file)
            list_of_genomes.append(n_name)

            with open(new_file, "w") as output_fh:
                for record in SeqIO.parse(file, "fasta"):
                    output_fh.write(">lcl|" + n_name + "|" + record.description + "\n")
                    output_fh.write(str(record.seq) + "\n")
        else:
            list_of_files.append(file)
            list_of_genomes.append(genome_name)

    return list_of_genomes, list_of_files


def get_genome_name(header):
    """
    Getting the name of the genome by hierarchy. This requires reading the first
    fasta header from the file. It also assumes a single genome per file.

    :param header: The header containing the record.
    :return genomeName: Name of the genome contained in the header.
    """

    re_patterns = (
        # Look for lcl followed by the possible genome name
        re.compile('(lcl\|[\w\-\.]+)'),

        # Look for contigs in the wwwwdddddd format
        re.compile('([A-Za-z]{4}\d{2})\d{6}'),

        # Look for a possible genome name at the beginning of the record ID
        re.compile('^(\w{8}\.\d)'),

        # Look for ref, gb, emb or dbj followed by the possible genome name
        re.compile('(ref\|\w{2}_\w{6}|gb\|\w{8}|emb\|\w{8}|dbj\|\w{8})'),

        # Look for gi followed by the possible genome name
        re.compile('(gi\|\d{8})'),


        # Look for name followed by space, then description
        re.compile('^([\w\-\.]+)\s+[\w\-\.]+')
    )

    # if nothing matches, use the full header as genome_name
    genome_name = header
    for rep in re_patterns:
        m = rep.search(header)

        if m:
            genome_name = m.group(1)
            break

    return str(genome_name)


def get_fasta_header_from_file(filename):
    """
    Gets the first fasta sequence from the file, and returns the fasta header.
    The files should have already been validated as fasta format.

    :param filename: the absolute path of the fasta file
    :return: header
    """

    for record in SeqIO.parse(filename, "fasta"):
        return record.description


def get_parsing_dict(ptype):
    """
    Given the parsed arguments from argparser, return a dictionary of functions.
    :param ptype: Type of parsing dict to return
    :return: {parser: function, predictor: function, data: data, type: type}
    """

    if ptype == 'serotype':
        # We will attach the JSON of known fasta headers and alleles to
        # a 'data' key in the parsing dictionary.
        json_handle = open(definitions.SEROTYPE_ALLELE_JSON, 'r')
        parsing_dict = {
            'parser':serotypePrediction.parse_serotype,
            'predictor':serotypePrediction.predict_serotype,
            'data':json.load(json_handle),
            'type':'serotype'}
        json_handle.close()
        return parsing_dict
    else:
        log.error("No parsing dictionary assigned for {0}".format(ptype))
        exit(1)

def assemble_reads(reads, reference):
    '''
    Return path to assembled reads.
    Assemble short reads by mapping to a reference genome.
    Default output is the same as reads file
        (basename+'iden.fasta' and basename+'pred.fasta').

    Args:
        reads (str): FASTQ/FQ format reads file
        reference (str): FASTA format reference file

    Returns:
        tuple(str, str): identifcation and prediction fasta file
    '''
    temp_dir = tempfile.gettempdir()
    output = os.path.join(
        temp_dir,
        os.path.splitext(os.path.basename(reads))[0]+'.fasta'
    )

    # run once if index reference does not exist
    index_path = \
        os.path.join(
            temp_dir,
            'bowtie_index',
            os.path.splitext(os.path.basename(reference))[0],
            'index'
        )
    index_dir = os.path.split(index_path)[0]
    if not os.path.isfile(index_path+'1.bt2'):
        log.info('Reference index does not exist. Creating new reference index at %s', index_dir)
        if not os.path.exists(index_dir):
            os.makedirs(index_dir)
        cmd1 = [
            'bowtie2-build',
            reference,
            index_path
        ]
        subprocess_util.run_subprocess(cmd1)

    cmd2 = [
        'bowtie2',
        '--score-min L,1,-0.5',
        '--np 5',
        '-x', index_path,
        '-U', reads,
        '-S', os.path.join(temp_dir, 'reads.sam')
    ]
    subprocess_util.run_subprocess(cmd2)

    cmd3 = [
        definitions.SAMTOOLS, 'view',
        '-F 4',
        '-q 1',
        '-bS', os.path.join(temp_dir, 'reads.sam'),
        '-o', os.path.join(temp_dir, 'reads.bam')
    ]
    subprocess_util.run_subprocess(cmd3)
    cmd4 = [
        definitions.SAMTOOLS, 'sort',
        os.path.join(temp_dir, 'reads.bam'),
        '-o', os.path.join(temp_dir, 'reads.sorted.bam'),
    ]
    subprocess_util.run_subprocess(cmd4)

    shell_cmd = [
        definitions.SAMTOOLS+' mpileup -uf', # mpileup
        reference,
        os.path.join(temp_dir, 'reads.sorted.bam'),
        '|',
        'bcftools call -c', # variant calling
        '|',
        'vcfutils.pl vcf2fq', # vcf to fq
        '|',
        'seqtk seq -A -', # fq to fasta
        '>',
        output
    ]
    subprocess_util.run_subprocess(' '.join(shell_cmd), is_shell=True)
    return split_mapped_output(output)

def get_num_hits(target, reference, args):
    '''
    Return number of matching hits when query the reference genome
        on the target genome
    
    Args:
        target(str): target genome
        reference(str): reference genome
        args (arguments): console arguments

    Returns:
        int: number of hits found
    '''
    num_hit = 0
    blast_db = blastFunctions.create_blast_db([target])
    result = blastFunctions.run_blast(
        definitions.ECOLI_MARKERS,
        blast_db,
        args
    )
    temp_dict = {}
    with open(result) as handler:
        log.debug("get_num_hits() output:")
        for line in handler:
            clean_line = line.strip()
            log.debug(clean_line)
            la = clean_line.split()
            temp_dict[la[0]] = True
    num_hit = len(temp_dict)
    log.debug("%s aligned to %d marker sequences", target, num_hit)
    return num_hit

def is_ecoli_genome(file, args):
    '''
    Return True if file is classified as ecoli by ecoli markers, otherwise False

    Args:
        file (str): path to valid fasta/fastq genome file
        args (arguments): console arguments

    Returns:
        bool: output
    '''
    num_hit = get_num_hits(file, definitions.ECOLI_MARKERS, args)
    if num_hit < 3:
        log.info("%s is not a valid e.coli genome file", file)
        return False
    log.info("%s is a valid e.coli genome file", file)
    return True

def get_num_of_fasta_entry(file):
    '''
    Return number of entries in a fasta file

    Args:
        file(str): path to fasta file

    Returns:
        int: number of entries
    '''
    count =0
    for _ in SeqIO.parse(file, 'fasta'):
        count += 1
    return count

def split_mapped_output(file):
    '''
    Split given fasta file into two file based on 'lcl' tags
        in the seq header
    Args:
        file(str): path to input fasta file

    Returns: 
        (str): path to ecoli identification fasta seq
        (str): path to serotype prediction fasta seq
    '''
    identif_file = os.path.splitext(file)[0]+'.iden.fasta'
    predict_file = os.path.splitext(file)[0]+'.pred.fasta'
    identif_seqs = []
    predict_seqs = []
    for record in SeqIO.parse(file, 'fasta'):
        if 'lcl' in record.description:
            identif_seqs.append(record)
        else:
            predict_seqs.append(record)
    with open(identif_file, "w") as output_handle:
        SeqIO.write(identif_seqs, output_handle, "fasta")
    with open(predict_file, "w") as output_handle:
        SeqIO.write(predict_seqs, output_handle, "fasta")
    return identif_file, predict_file