import sys
import tempfile
import os
import logging
from ectyper import ectyper


TEST_ROOT = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(TEST_ROOT, '..'))
MASHINFILE = os.path.join(ROOT_DIR, 'ectyper/Data/refseq.genomes.k21s1000.msh')


def set_input(input,
              percent_iden=None,
              verify=True,
              output=tempfile.mkdtemp(),
              cores=1,
              print_sequence=False,
              refseqmash = False):
    """
    Create the sys.argv[] without need for commandline input.
    :param input: Input file given by testing function
    :param percent_iden: Percent identity for comparison
    :param output: Location of output
    :return: None
    """
    args = ['-i', input,
            '-c', str(cores),
            ]

    if percent_iden:
        args += ['-d', str(percent_iden)]
    if verify:
        args += ['--verify']
    if output:
        args += ['-o', output]
    if print_sequence:
        args += ['--sequence']
    if refseqmash:
        args += ['-r', MASHINFILE]

    sys.argv[1:] = args



def test_failed_species_identification(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT, 'Data/GCF_001672015.1.fna')
    set_input(input=file, verify=True)
    ectyper.run_program()
    assert "GCF_001672015.1\tEscherichia coli\t-\tH8\t-:H8\tWARNING (-:H TYPING)" in caplog.text

def test_failed_species_identification_nospeciesverify(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT, 'Data/GCF_001672015.1.fna')
    set_input(input=file, verify=False)
    ectyper.run_program()
    assert "GCF_001672015.1\t-\t-\tH8\t-:H8\t-" in caplog.text


def test_non_existing_accession_in_meta(caplog):
    """
    GCA_900059685.2 - Streptococcus pneumoniae - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/059/685/GCA_900059685.2_12291_5_44
    is not present in assembly_summary_refseq.txt and is a perfect candidate to try to test
    species identification function
    :param caplog:
    :return:
    """
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT, 'Data/GCA_900059685.2.fna')
    set_input(input=file, verify=False)
    ectyper.run_program()
    assert "No O and H antigen determinant E.coli genes were found" in caplog.text

def test_speciesID_non_existing():
    """
    GCF_001672015.1.fna is not present in `assembly_summary_refseq.txt`. Let's see how get_species would recover from this
    using top 10 hits
    :return:
    """
    fastafile=os.path.join(TEST_ROOT, 'Data/GCF_001672015.1.fna')
    set_input(input=fastafile, verify=False, refseqmash=True)
    args = ectyper.commandLineOptions.parse_command_line()
    species = ectyper.speciesIdentification.get_species(file=fastafile, args=args)
    print("Identified species {}".format(species))
    assert "Escherichia coli" == species