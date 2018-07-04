#!/usr/bin/env python

'''
Genome Utilities
'''

import logging
import os
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
    if file_or_directory:
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

    if not files_list:
        LOG.critical("No files were found for the ectyper run")
        exit("No files were found")

    sorted_files = sorted(files_list)
    LOG.debug(sorted_files)
    return sorted_files


def get_valid_format(file):
    """
    Check using SeqIO if files are valid fasta/fastq format, returns the format.

    Args:
        file (str): path of file

    Returns:
        fm (str or None): the file format if 'fasta', 'fastq', 'other'
    """
    for fm in ['fastq', 'fasta']:
        try:
            with open(file, "r") as handle:
                data = SeqIO.parse(handle, fm)
                if any(data):
                    if is_tarfile(file):
                        LOG.warning("Compressed file is not supported: {}".format(file))
                        return 'other'
                    return fm
        except FileNotFoundError:
            LOG.warning("{0} is not found.".format(file))
            return 'other'
        except UnicodeDecodeError:
            LOG.warning("{0} is not a valid file.".format(file))
            return 'other'
        except ValueError:
            LOG.debug("{0} is not a {1} file.".format(file, fm))

    LOG.warning("{0} is not a fasta/fastq file".format(file))
    return 'other'


def get_genome_names_from_files(files, temp_dir):
    """
    For each file:
    Uses the name of the file for the genome name, creates a temporary file using
    >lcl|filename as the name in the fasta header.
    :param files: All the fasta files for analyses
    :param temp_dir: The ectyper temp directory
    :return: List of files with fasta headers modified for filename
    """

    modified_genomes = []
    for file in files:
        # get only the name of the file for use in the fasta header
        file_base_name = os.path.basename(file)
        file_path_name = os.path.splitext(file_base_name)[0]
        n_name = file_path_name.replace(' ', '_')

        # create a new file for the updated fasta headers
        new_file = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False).name

        # add the new name to the list of files and genomes
        modified_genomes.append(new_file)

        with open(new_file, "w") as outfile:
            with open(file) as infile:
                for record in SeqIO.parse(infile, "fasta"):
                    outfile.write(">lcl|" + n_name + "|" + record.description + "\n")
                    outfile.write(str(record.seq) + "\n")

    LOG.debug(("Modified genomes: {}".format(modified_genomes)))
    return modified_genomes


def create_bowtie_base(temp_dir, reference):
    """
    Create the bowtie reference, based on the combined E. coli markers and
    O- and H- type alleles

    :param temp_dir: Program wide temporary directory
    :param reference: The fasta file to use as a reference for the base
    :return: The full path to bowtie_base
    """

    LOG.info("Creating bowtie2 base index")

    bowtie_base = os.path.join(temp_dir, 'bowtie_reference')
    LOG.info("Creating the bowtie2 index at {}".format(bowtie_base))
    bowtie_build = [
        'bowtie2-build',
        '-f',
        reference,
        bowtie_base
    ]
    subprocess_util.run_subprocess(bowtie_build)
    return bowtie_base


def assemble_reads(reads, bowtie_base, combined_fasta, temp_dir):
    """
    Assembles fastq reads to the specified reference file.
    :param reads: The fastq file to assemble
    :param bowtie_base: The full-path to the bowtie reference created for this run
    :param combined_fasta: Combined E. coli markers and O- and H- alleles
    :param temp_dir: The ectyper temporary directory
    :return: The full path to the assembled fasta file
    """

    output_fasta = os.path.join(
        temp_dir,
        os.path.splitext(os.path.basename(reads))[0] + '.fasta'
    )

    # Run bowtie2
    sam_reads = os.path.join(temp_dir, 'reads.sam')
    bowtie_run = [
        'bowtie2',
        '--score-min L,1,-0.5',
        '--np 5',
        '-x', bowtie_base,
        '-U', reads,
        '-S', sam_reads
    ]
    subprocess_util.run_subprocess(bowtie_run)

    # Convert reads from sam to bam
    bam_reads = os.path.join(temp_dir, 'reads.bam')
    sam_convert = [
        'samtools', 'view',
        '-S',
        '-F 4',
        '-q 1',
        '-b',
        '-o', bam_reads,
        sam_reads,
    ]
    subprocess_util.run_subprocess(sam_convert)

    # Sort the reads
    sorted_bam_reads = os.path.join(temp_dir, 'reads.sorted.bam')
    sam_sort = [
        'samtools', 'sort',
        bam_reads,
        '-o', sorted_bam_reads
    ]
    subprocess_util.run_subprocess(sam_sort)

    # Create fasta from the reads
    mpileup = [
        'bcftools', 'mpileup',
        '-f', combined_fasta,
        sorted_bam_reads]
    mpileup_output = subprocess_util.run_subprocess(mpileup)

    variant_calling = [
        'bcftools',
        'call',
        '-c'
    ]
    variant_calling_output = subprocess_util.run_subprocess(variant_calling, mpileup_output.stdout)

    to_fastq = [
        'vcfutils.pl',
        'vcf2fq'
    ]
    to_fastq_output = subprocess_util.run_subprocess(to_fastq, variant_calling_output.stdout)

    to_fasta = [
        'seqtk',
        'seq',
        '-A'
    ]
    to_fasta_output = subprocess_util.run_subprocess(to_fasta, to_fastq_output.stdout)

    # Write the final fasta file as bytes rather than str
    LOG.info("Creating fasta file {}".format(output_fasta))
    with open(output_fasta, 'wb') as ofh:
        ofh.write(to_fasta_output.stdout)

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
    other_files = []

    for file in raw_files:
        file_format = get_valid_format(file)
        if file_format == 'fasta':
            fasta_files.append(file)
        elif file_format == 'fastq':
            fastq_files.append(file)
        else:
            other_files.append(file)

    LOG.debug('raw fasta files: {}'.format(fasta_files))
    LOG.debug('raw fastq files: {}'.format(fastq_files))
    LOG.debug("other non- fasta/fastq files: {}".format(other_files))

    return({'fasta':fasta_files,
            'fastq':fastq_files,
            'other':other_files})


def assembleFastq(raw_files_dict, temp_dir, combined_fasta, bowtie_base):
    """
    For any fastq files, map and assemble the serotyping genes, and optionally
    the E. coli specific genes.
    :param raw_files_dict: Dictionary of ['fasta'] and ['fastq'] files
    :param temp_dir: Temporary files created for ectyper
    :param combined_fasta: Combined E. coli markers and O- and H- alleles
    :param bowtie_base: The bowtie base index of O- and H- alleles and E. coli markers
    :return: list of all fasta files, including the assembled fastq
    """

    all_fasta_files = raw_files_dict['fasta']
    for fastq_file in raw_files_dict['fastq']:
        fasta_file = assemble_reads(fastq_file, bowtie_base, combined_fasta, temp_dir)
        all_fasta_files.append(fasta_file)

    return all_fasta_files


def create_combined_alleles_and_markers_file(alleles_fasta, temp_dir):
    """
    Create a combined E. coli specific marker and serotyping allele file
    Do this every program run to prevent the need for keeping a separate combined file in sync.

    :param alleles_fasta: The O- and H- alleles file
    :param temp_dir: Temporary directory for program run
    :return: The combined alleles and markers file
    """

    combined_file = os.path.join(temp_dir, 'combined_ident_serotype.fasta')
    LOG.info("Creating combined serotype and identification fasta file")

    with open(combined_file, 'w') as ofh:
        with open(definitions.ECOLI_MARKERS, 'r') as mfh:
            ofh.write(mfh.read())

        with open(alleles_fasta, 'r') as sfh:
            ofh.write(sfh.read())

    return combined_file
