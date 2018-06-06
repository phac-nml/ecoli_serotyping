import sys
import pytest
import os
# sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
# from ectyper import ectyper
from ectyper import ectyper

def set_input(input, percent_iden=None, output=None):
    """
    Create the sys.argv[] without need for commandline input.
    :param input: Input file given by testing function
    :param percent_iden: Percent identity for comparison
    :param output: Location of output
    :return: None
    """
    args = ['-i', input]
    if percent_iden:
        args += ['-d', str(percent_iden)]
    if output:
        args += ['-o', output]
    sys.argv[1:] = args


def test_integration_invalid_file():
    """
    Giving a non-fasta file in fasta-file name.
    :return: None
    """
    file = 'test/Data/test_dir/badfasta.fasta'
    set_input(file)
    with pytest.raises(SystemExit) as se:
        ectyper.run_program()
    assert se.type == SystemExit
    assert se.value.code == "No valid genomes"


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
    assert se.value.code == "No valid genomes"

