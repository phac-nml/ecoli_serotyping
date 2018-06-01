#!/usr/bin/env python

'''
Genome Utilities
'''

import logging
import os
import re
import tempfile
from tarfile import is_tarfile

from Bio import SeqIO
from ectyper import definitions, subprocess_util

LOG = logging.getLogger(__name__)


def get_files_as_list(file_or_directory):
    """
    Creates a list of files from either the given file, or all files within the
    directory specified (where each file name is its absolute path).

    Args:
        file_or_directory (str): file or directory name given on commandline

    Returns:
        files_list (list(str)): List of all the files found.

    """

    files_list = []
    if file_or_directory == '':
        return files_list

    if os.path.isdir(file_or_directory):
        LOG.info("Gathering genomes from directory " + file_or_directory)

        # Create a list containing the file names
        for root, dirs, files in os.walk(file_or_directory):
            for filename in files:
                files_list.append(os.path.join(root, filename))
    # check if input is concatenated file locations
    elif ',' in file_or_directory:
        LOG.info("Using genomes in the input list")
        for filename in file_or_directory.split(','):
            files_list.append(os.path.abspath(filename))
    else:
        LOG.info("Using genomes in file " + file_or_directory)
        files_list.append(os.path.abspath(file_or_directory))

    return sorted(files_list)


def get_valid_format(file):
    """
    Check using SeqIO if files are valid fasta/fastq format, returns the format.

    Args:
        file (str): path of file

    Returns:
        fm (str or None): the file format if 'fasta' or 'fastq', otherwise None
    """
    for fm in ['fastq', 'fasta']:
        try:
            with open(file, "r") as handle:
                data = SeqIO.parse(handle, fm)
                if any(data):
                    if is_tarfile(file):
                        LOG.warning("Compressed file is not supported: {}".format(file))
                        return None
                    return fm
        except FileNotFoundError as err:
            LOG.warning("{0} is not found".format(file))
            return None
        except UnicodeDecodeError as err:
            LOG.warning("{0} is not a valid file".format(file))
            return None
        except:
            LOG.warning("{0} is an unexpected file".format(file))
            return None
    LOG.warning("{0} is not a fasta/fastq file".format(file))
    return None


def get_genome_names_from_files(files, temp_dir):
    """
    For each file:
    Takes the first header from a fasta file and sends it to the get_genome_name
    function. Returns the name of the genome. If the name of the file is to be
    used as the genome name, creates a temporary file using >lcl|filename as the
    first part of the header.

    Args:
        files (list): The list of files to get the genome names for
        temp_dir (str): A tempdir where the copied files will be stored

    Returns:
        tuple(list(str), list(str)): first list is genome names, second list is file names
    """

    list_of_genomes = []
    list_of_files = []
    for file in files:
        # Always to use the filename as the genome name.
        # This means we also need to create a new file adding
        # the filename to each of the headers, so that downstream applications
        # (eg. BLAST) can be used with the filename as genome name.

        # get only the name of the file for use in the fasta header
        file_base_name = os.path.basename(file)
        file_path_name = os.path.splitext(file_base_name)[0]
        n_name = file_path_name.replace(' ', '_')

        # create a new file for the updated fasta headers
        new_file = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False).name

        # add the new name to the list of files and genomes
        list_of_files.append(new_file)
        list_of_genomes.append(n_name)

        with open(new_file, "w") as outfile:
            with open(file) as infile:
                for record in SeqIO.parse(infile, "fasta"):
                    outfile.write(">lcl|" + n_name + "|" + record.description + "\n")
                    outfile.write(str(record.seq) + "\n")

    return list_of_genomes, list_of_files


def get_genome_name(header):
    """
    Getting the name of the genome by hierarchy. This requires reading the first
    fasta header from the file. It also assumes a single genome per file.

    Args:
        header (str): The header containing the record.

    Returns:
        genomeName (str): Name of the genome contained in the header.
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

def assemble_reads(reads, reference, temp_dir):
    '''
    Return path to assembled reads.
    Assemble short reads by mapping to a reference genome.
    Default output is the same as reads file
        (basename+'iden.fasta' and basename+'pred.fasta').

    Args:
        reads (str): FASTQ/FQ format reads file
        reference (str): FASTA format reference file
        temp_dir (str): temp_dir for storing assembled files

    Returns:
        tuple(str, str): identifcation and prediction fasta file
    '''
    output = os.path.join(
        temp_dir,
        os.path.splitext(os.path.basename(reads))[0]+'.fasta'
    )

    # run once if index reference does not exist
    index_path = \
        os.path.join(
            definitions.DATA_DIR,
            'bowtie_index',
            os.path.splitext(os.path.basename(reference))[0]
        )
    index_dir = os.path.split(index_path)[0]
    if not os.path.isdir(index_dir):
        LOG.info('Reference index does not exist. Creating new reference index at {}'.format(index_dir))
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

def split_mapped_output(file):
    '''
    Splits given fasta files into two file based on 'lcl' tags
        in the seq header

    Args:
        file (str): path to input fasta file

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
