import sys
import tempfile
import os
import logging
from ectyper import ectyper


TEST_ROOT = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(TEST_ROOT, '..'))
MASHINFILE = os.path.join(ROOT_DIR, 'ectyper/Data/EnteroRef_GTDBSketch_20231003_V2.msh')


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
    assert "GCF_001672015.1\tEscherichia coli" in caplog.text
    assert "-\tH8\t-:H8\tWARNING (-:H TYPING)" in caplog.text

def test_failed_species_identification_nospeciesverify(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT, 'Data/GCF_001672015.1.fna')
    set_input(input=file, verify=False)
    ectyper.run_program()
    assert "GCF_001672015.1\tEscherichia coli" in caplog.text
    assert "-\tH8\t-:H8\t-" in caplog.text


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
    GCF_001672015.1.fna species ID identification making sure it is identified as E.coli
    :return:
    """
    fastafile=os.path.join(TEST_ROOT, 'Data/GCF_001672015.1.fna')
    set_input(input=fastafile, verify=False, refseqmash=True)
    args = ectyper.commandLineOptions.parse_command_line()
    result = ectyper.speciesIdentification.verify_ecoli_and_inputs(fasta_fastq_files_dict={fastafile:None}, 
                                                                   ofiles={}, filesnotfound={}, args=args)
    assert result[0]['GCF_001672015.1']['species'] == 'Escherichia coli'

def test_speciesID_Shigella(caplog):
    fastafile=os.path.join(TEST_ROOT, 'Data/ESC_BA0255AA_AS_Shigella_sonnei.fasta')
    caplog.set_level(logging.DEBUG)
    tmpdir = tempfile.mkdtemp()
    set_input(input=fastafile, verify=False, output=tmpdir)
    expected_species = "Shigella sonnei"
    
    ectyper.run_program()
    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
         rows = outfp.readlines()
    species_col_number = [idx for idx, i in enumerate(rows[0].split('\t')) if i == "Species"]
    assert len(species_col_number) == 1, "Could not find the species column"     
    secondrow=rows[1:][0].split('\t') #check only second row
    assert secondrow[species_col_number[0]] == expected_species, f"Could not find species {expected_species}"

def test_Ealbertii_1(caplog): #error
    logging.info("Starting 1 of 3 test on EnteroBase on sample ESC_HA8355AA_AS: Escherichia albertii O65:H5")
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/ESC_HA8355AA_AS_Ealberii_O65H5.fasta')
    tmpdir = tempfile.mkdtemp()
    set_input(input=file, cores=4, print_sequence=True,  verify=True, output=tmpdir)
    ectyper.run_program()

    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
         rows = outfp.readlines()
    secondrow=rows[1:][0] #remove header line
    assert "Escherichia albertii" in secondrow
    assert "WARNING (WRONG SPECIES)" in secondrow

def test_Ealbertii_2(caplog): #error
    logging.info("Starting 2 of 3 test on EnteroBase on sample on ESC_HA8509AA_AS: Escherichia albertii O5:H5")
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/ESC_HA8509AA_AS_EalbertiiO5H5.fasta')
    tmpdir = tempfile.mkdtemp()
    set_input(input=file, cores=4, print_sequence=True, verify=True, output=tmpdir)
    ectyper.run_program()

    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
         rows = outfp.readlines()
    secondrow=rows[1:][0] #check only second row

    assert "Escherichia albertii" in secondrow
    assert "WARNING (WRONG SPECIES)" in secondrow

def test_Ealbertii_3(caplog):
    logging.info("Starting 3 of 3 test Escherichia albertii O49:NM") #can not type O49 due to poor sequence quality of uncertainty of wet-lab O49 typing
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/Ealbertii_O49NM.fasta')

    tmpdir = tempfile.mkdtemp()
    set_input(input=file, cores=4, print_sequence=True,  verify=True, output=tmpdir)
    ectyper.run_program()

    with open(os.path.join(tmpdir ,"output.tsv")) as outfp:
         rows = outfp.readlines()
    secondrow=rows[1:][0] #check only second row
    assert "Escherichia albertii" in secondrow
    assert "WARNING (WRONG SPECIES)" in secondrow