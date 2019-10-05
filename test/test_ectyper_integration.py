import sys
import pytest
import tempfile
import os
import logging
from ectyper import ectyper

TEST_ROOT = os.path.dirname(__file__)


def set_input(input,
              percent_iden=None,
              output=tempfile.mkdtemp(),
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
            '--verify',
            '-r', os.path.join(TEST_ROOT, 'Data/test_sketch.msh'),
            '-c', str(cores)]

    if percent_iden:
        args += ['-d', str(percent_iden)]
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
    set_input(file)
    ectyper.run_program()
    assert "Non fasta / fastq file" in caplog.text


def test_integration_no_file():
    """
    Giving no input to the program.
    :return: None
    """
    file = ''
    set_input(file)
    with pytest.raises(SystemExit) as se:
        ectyper.run_program()
    assert se.type == SystemExit
    assert se.value.code == "No files were found"


def test_integration_valid_file(caplog):
    """
    Ensure a valid E. coli fasta passes
    :return: None
    """
    file = os.path.join(TEST_ROOT, 'Data/Escherichia.fna')
    set_input(file)
    ectyper.run_program()
    assert "Escherichia\tEscherichia coli\tO103\tH2\tO103:H2\tPASS\tHIGH" in caplog.text


def test_integration_yersinia(caplog):
    """
    Ensure a non-E. coli gets categorized as such
    :return: None
    """
    file = os.path.join(TEST_ROOT, 'Data/Yersinia.fasta')
    set_input(file)
    ectyper.run_program()
    assert "Yersinia pestis" in caplog.text


def test_valid_fastq_file(caplog):
    """
    Given a valid fastq file, get the correct results.
    Use a temp dir for the test output
    :return: None
    """
    file = os.path.join(TEST_ROOT, 'Data/Escherichia.fastq')
    set_input(file)
    ectyper.run_program()
    assert "Escherichia\tEscherichia coli\tO22\tH8" in caplog.text


def test_multiple_directories(caplog):
    """
    Check a number of small files, some good, some bad,
    within a nested directory structure.

    :param caplog: Capture logging output for pytest
    :return: None
    """
    the_dir = os.path.join(TEST_ROOT, 'Data/test_dir')
    set_input(the_dir, cores=4, print_sequence=True)
    ectyper.run_program()
    assert "sample2\tEscherichia coli\tO148\tH44" in caplog.text
    assert "sample3\tEscherichia coli\tO148\tH44" in caplog.text
    assert "sample4\tEscherichia coli\tO148\tH44" in caplog.text
    assert "badfasta\t-\t-\t-\t-:-\t-\t-\t-\t-\tNon fasta / fastq file" in caplog.text
    assert "sample.fasta\t-\t-\t-\t-:-\t-\t-\t-\t-\tNon fasta / fastq file" in caplog.text
    assert "sampletar\t-\t-\t-\t-:-\t-\t-\t-\t-\tNon fasta / fastq file" in caplog.text
    assert "test_junk\t-\t-\t-\t-:-\t-\t-\t-\t-\tNon fasta / fastq file" in caplog.text

