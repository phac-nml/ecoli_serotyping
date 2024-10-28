import sys
import pytest
import tempfile
import os
from ectyper import ectyper, definitions
import subprocess
import pandas
import logging
import re

TEST_ROOT = os.path.dirname(__file__)
LOG=logging.getLogger("TEST")
LOG.setLevel(logging.INFO)


def set_input(input,
              percent_iden=None,
              verify = True,
              output=tempfile.mkdtemp(),
              maxdirdepth=0,
              cores=1,
              print_sequence=False):
    """
    Create the sys.argv[] without need for commandline input.
    :param input: Input file given by testing function
    :param percent_iden: Percent identity for comparison
    :param output: Location of output
    :return: None
    """
    args = ['-i', input,
            '--maxdirdepth', str(maxdirdepth),
            '-c', str(cores)]

    if percent_iden:
        args += ['-d', str(percent_iden)]
    if verify:
        args += ['--verify']
    if output:
        args += ['-o', output]
    if print_sequence:
        args += ['--sequence']
    sys.argv[1:] = args


def test_integration_invalid_file(caplog):
    """
    Giving a non-fasta file in fasta-file name.
    :return: None
    """
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT, 'Data/test_dir/badfasta.fasta')
    set_input(input=file)
    ectyper.run_program()
    assert "Non fasta / fastq file" in caplog.text


def test_integration_no_file():
    """
    Giving no input to the program.
    :return: None
    """
    file = ''
    set_input(input=file)
    with pytest.raises(FileNotFoundError) as se:
        ectyper.run_program()
    assert se.type == FileNotFoundError
    assert str(se.value) == "No files were found to run on"


def test_integration_valid_file(caplog):
    """
    Ensure a valid E. coli fasta passes
    :return: None
    """
    file = os.path.join(TEST_ROOT, 'Data/Escherichia.fna')
    set_input(file)
    ectyper.run_program()
    print(caplog.text)
    assert "PASS (REPORTABLE)" in caplog.text
    assert "O103:H2" in caplog.text
    assert "Escherichia coli" in caplog.text


def test_integration_yersinia(caplog):
    """
    Ensure a non-E. coli gets categorized as such
    :return: None
    """
    file = os.path.join(TEST_ROOT, 'Data/Yersinia.fasta')
    set_input(file)
    ectyper.run_program()
    assert "Yersinia pestis" in caplog.text
    assert "WARNING (WRONG SPECIES)" in caplog.text

def test_integration_validfasta_noverify(caplog):
    """
    Tests for fasta files without E.coli species verify function (--verify) do not fail as per issue #76 (https://github.com/phac-nml/ecoli_serotyping/issues/76)
    :return: None
    """
    file = os.path.join(TEST_ROOT, 'Data/Escherichia.fna')
    set_input(file, verify=False)
    ectyper.run_program()
    print(caplog.text)
    assert "O103\tH2\tO103:H2" in caplog.text
    assert "Escherichia coli" in caplog.text

def test_valid_fastq_file(caplog):
    """
    Given a valid fastq file, get the correct results.
    Use a temp dir for the test output
    :return: None
    """
    file = os.path.join(TEST_ROOT, 'Data/Escherichia.fastq')
    set_input(file, verify=False)
    ectyper.run_program()
    assert "O22:H8" in caplog.text

def test_valid_fastq_file_with_verify(caplog):
    """
    Given a valid fastq file with low genome coverage, test species verification fail
    Use a temp dir for the test output
    :return: None
    """
    file = os.path.join(TEST_ROOT, 'Data/Escherichia.fastq')
    set_input(file, verify=True)
    ectyper.run_program()
    assert "Escherichia coli" in caplog.text

def test_multiple_directories(caplog):
    """
    Check a number of small files, some good, some bad,
    within a nested directory structure.

    :param caplog: Capture logging output for pytest
    :return: None
    """
    the_dir = os.path.join(TEST_ROOT, 'Data/test_dir')
    set_input(the_dir, cores=4, verify=True, maxdirdepth=1, print_sequence=True)
    ectyper.run_program()
    assert any([True if re.match(r".+sample2.+WARNING\s+\(WRONG\s+SPECIES\)", line) else False for line in caplog.text.splitlines()]), "Issue with sample 2"
    assert any([True if re.match(r".+sample3.+WARNING\s+\(WRONG\s+SPECIES\)", line) else False for line in caplog.text.splitlines()]), "Issue with sample 3"
    assert any([True if re.match(r".+sample4.+WARNING\s+\(WRONG\s+SPECIES\)", line) else False for line in caplog.text.splitlines()]), "Issue with sample 4" 
    assert any([True if re.match(r".+badfasta.+WARNING\s+\(WRONG\s+SPECIES\).+Non fasta / fastq file", line) else False for line in caplog.text.splitlines()]), "Issue with badfasta"
    assert any([True if re.match(r".+sample.fasta.+WARNING\s+\(WRONG\s+SPECIES\).+Non fasta / fastq file", line) else False for line in caplog.text.splitlines()]), "Issue with sample.fasta"
    assert any([True if re.match(r".+sampletar.+WARNING\s+\(WRONG\s+SPECIES\).+Non fasta / fastq file", line) else False for line in caplog.text.splitlines()]), "Issue with sampletar"
    assert any([True if re.match(r".+test_junk.+WARNING\s+\(WRONG\s+SPECIES\).+Non fasta / fastq file", line) else False for line in caplog.text.splitlines()]), "Issue with test_junk"
    assert any([True if re.match(r".+GCA_000181775\.1_ASM18177v1_genomic\s+Escherichia\s+coli.+O157\s+H7.+REPORTABLE", line) else False for line in caplog.text.splitlines()]), "Issue with GCF_000181775.1_ASM18177v1"


def test_mash_sketch_and_assembly_metadata(tmpdir):
    """
    Test if all accessions in mash sketch are a complete subset of the assembly stats superset.
    Ensure that all accession numbers are represented in the meta data assembly stats
    """
    ectyper.speciesIdentification.get_species_mash(definitions.SPECIES_ID_SKETCH)
    MASHSTATSMETAFILE=os.path.join(tmpdir+"/species_id_mash_meta.txt")
    MASHINFILE = definitions.SPECIES_ID_SKETCH

    cmd = ["mash info -t " +  MASHINFILE   + " > " + MASHSTATSMETAFILE]
    LOG.info("Mash info text file written to {}".format(MASHSTATSMETAFILE))
    subprocess.run(cmd, shell=True)
    mashsketchdatadf = pandas.read_csv(MASHSTATSMETAFILE,sep="\t")
    assert mashsketchdatadf.columns.to_list() == ['#Hashes', 'Length', 'ID', 'Comment']
    assert ['GCA' in item for item in mashsketchdatadf["ID"]], "GCA_ identifiers not found"
    assert ['s__Escherichia coli' in item for item in mashsketchdatadf["Comment"]], "GCA_ identifiers not found"





