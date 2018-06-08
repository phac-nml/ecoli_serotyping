#!/usr/bin/env python

'''
Genome Utilities
'''

import logging
import os
import sys
import tempfile
from tarfile import is_tarfile
from Bio import SeqIO
from ectyper import definitions, subprocess_util, speciesIdentification
from urllib.request import urlretrieve
import shutil
from collections import defaultdict

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

    sorted_files = sorted(files_list)
    LOG.debug(sorted_files)
    return sorted_files


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
            except FileNotFoundError:
                LOG.warning("{0} is not found.".format(file))
                return None
            except UnicodeDecodeError:
                LOG.warning("{0} is not a valid file.".format(file))
                return None
            except ValueError:
                LOG.debug("{0} is not a {1} file.".format(file, fm))
                return None

    LOG.warning("{0} is not a fasta/fastq file".format(file))
    return None



def get_genome_names_from_files(files, temp_dir):
    """
    For each file:
    Uses the name of the file for the genome name, creates a temporary file using
    >lcl|filename as the name in the fasta header.

    Args:
        files (list): The list of files to get the genome names for
        temp_dir (str): A tempdir where the copied files will be stored

    Returns:
        tuple(list(str), list(str)): first list is genome names, second list is file names
    """

    list_of_genomes = []
    list_of_files = []
    for file in files:
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


def assemble_reads(reads, reference, temp_dir):
    """
    Assembles fastq reads to the specified reference file.
    :param reads: The fastq file to assemble
    :param reference: The fasta reference file
    :param temp_dir: The ectyper temporary directory
    :return: The full path to the assembled fasta file
    """

    output_fasta = os.path.join(
        temp_dir,
        os.path.splitext(os.path.basename(reads))[0],
        '.fasta'
    )

    # Create the bowtie2 reference index
    bowtie_base = os.path.join(temp_dir, 'bowtie_reference')
    bowtie_build = [
        'bowtie2-build',
        reference,
        bowtie_base
    ]
    subprocess_util.run_subprocess(bowtie_build)

    # Run bowtie2
    bowtie_run = [
        'bowtie2',
        '--score-min L,1,-0.5',
        '--np 5',
        '-x', bowtie_base,
        '-U', reads,
        '-S', os.path.join(temp_dir, 'reads.sam')
    ]
    subprocess_util.run_subprocess(bowtie_run)

    # Convert reads from sam to bam
    sam_convert = [
        definitions.SAMTOOLS, 'view',
        '-F 4',
        '-q 1',
        '-bS', os.path.join(temp_dir, 'reads.sam'),
        '-o', os.path.join(temp_dir, 'reads.bam')
    ]
    subprocess_util.run_subprocess(sam_convert)

    # Sort the reads
    sam_sort = [
        definitions.SAMTOOLS, 'sort',
        os.path.join(temp_dir, 'reads.bam'),
        '-o', os.path.join(temp_dir, 'reads.sorted.bam'),
    ]
    subprocess_util.run_subprocess(sam_sort)

    # Create fasta from the reads
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
        output_fasta
    ]
    subprocess_util.run_subprocess(' '.join(shell_cmd), is_shell=True)

    return output_fasta


def get_raw_files(raw_files):
    """Take all the raw files, and filter not fasta / fastq

    Args:
        raw_files(str): list of files from user input

    Returns:
        A dictionary collection of fasta and fastq files
        example:
        {'raw_fasta_files':[],
         'raw_fastq_files':[]}
    """
    fasta_files = []
    fastq_files = []

    for file in raw_files:
        file_format = get_valid_format(file)
        if file_format == 'fasta':
            fasta_files.append(file)
        elif file_format == 'fastq':
            fastq_files.append(file)

    LOG.debug('raw fasta files: {}'.format(fasta_files))
    LOG.debug('raw fastq files: {}'.format(fastq_files))

    return({'fasta':fasta_files, 'fastq':fastq_files})


def download_refseq():
    '''Download refseq file with progress bar
    '''

    def reporthook(blocknum, blocksize, totalsize):
        '''
        https://stackoverflow.com/questions/15644964/python-progress-bar-and-downloads
        '''
        readsofar = blocknum * blocksize
        if totalsize > 0:
            s = "\r {:5.1%} {:{}d} / {:d}".format(
                readsofar/totalsize, readsofar,
                len(str(totalsize)),
                totalsize
            )
            sys.stderr.write(s)
            if readsofar >= totalsize: # near the end
                sys.stderr.write("\n")
        else: # total size is unknown
            sys.stderr.write("read {}\n".format(readsofar))

    if not os.path.isfile(definitions.REFSEQ_SKETCH):
        refseq_url = 'https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh'
        LOG.info("No refseq found. Downloading reference file for species identification...")
        urlretrieve(refseq_url, definitions.REFSEQ_SKETCH, reporthook)
        LOG.info("Download complete.")


def assembleFastq(raw_files_dict, temp_dir):
    """
    For any fastq files, map and assemble the serotyping genes, and optionally
    the E. coli specific genes.
    :param raw_files_dict: Dictionary of ['fasta'] and ['fastq'] files
    :param temp_dir: Temporary files created for ectyper
    :return: list of all fasta files, including the assembled fastq
    """

    # Create a combined E. coli specific marker and serotyping allele file
    # Do this every program run to prevent the need for keeping a separate
    # combined file in sync.
    combined_file = os.path.join(temp_dir, 'combined_ident_serotype.fasta')
    shutil.copy(definitions.ECOLI_MARKERS, combined_file)
    shutil.copy(definitions.SEROTYPE_FILE, combined_file)

    all_fasta_files = raw_files_dict['fasta']
    for fastq_file in raw_files_dict['fastq']:
        fasta_file = assemble_reads(fastq_file, combined_file, temp_dir)
        all_fasta_files.append(fastq_file)

    return all_fasta_files




def filter_file_by_species(genome_file, genome_format, temp_dir, verify=False, species=False):
    """
    Assemble fastq sequences to fasta for the E. coli specific markers
    if verify is enabled. If identified as non-E. coli, identify the probable species
    using MASH, if enabled.

    Args:
        genome_file: input genome file
        genome_format(str): fasta or fastq
        temp_dir: temporary directory
        verify(bool):
            whether to perform E. coli verification
        species(bool):
            whether to perform species identification for non-E. coli genomes
    Returns:
        The filtered and assembled genome files in fasta format
    """

    filtered_file = None
    if genome_format == 'fastq':
        iden_file, pred_file = \
            assemble_reads(genome_file, definitions.ECOLI_MARKERS, temp_dir)
        # If no alignment result, the file is definitely not E. coli
        if get_valid_format(iden_file) is None:
            LOG.warning(
                "{} is filtered out because no identification alignment found".format(genome_file))
            return filtered_file
        if not (verify or species) or speciesIdentification.is_ecoli_genome(
                iden_file, genome_file, mash=species):
            # final check before adding the alignment for prediction
            if get_valid_format(iden_file) != 'fasta':
                LOG.warning(
                    "{0} is filtered out because no prediction alignment found".format(genome_file))
                return filtered_file
            filtered_file = pred_file

    if genome_format == 'fasta':
        if not (verify or species) \
        or speciesIdentification.is_ecoli_genome(genome_file, mash=species):
            filtered_file = genome_file
    return filtered_file
