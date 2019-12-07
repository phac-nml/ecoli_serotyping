import sys
import pytest
import tempfile
import os
from ectyper import ectyper
import subprocess
import pandas
import logging
import re

TEST_ROOT = os.path.dirname(__file__)
logging.basicConfig(level=logging.INFO)
LOG=logging.getLogger("TEST")


def set_input(input,
              percent_iden=None,
              verify = True,
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
            #'-r', os.path.join(TEST_ROOT, 'Data/test_sketch.msh'),
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

def test_integration_validfasta_noverify(caplog):
    """
    Tests for fasta files without E.coli species verify function (--verify) do not fail as per issue #76 (https://github.com/phac-nml/ecoli_serotyping/issues/76)
    :return: None
    """
    file = os.path.join(TEST_ROOT, 'Data/Escherichia.fna')
    set_input(file, verify=False)
    ectyper.run_program()
    assert "Escherichia\tEscherichia coli O103:H2 str. 12009\tO103\tH2\tO103:H2" in caplog.text


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

def test_mash_sketch_and_assembly_metadata():
    """
    Test if all accessions in mash sketch are a complete subset of the assembly stats superset.
    Ensure that all accession numbers are represented in the meta data assembly stats
    """
    ectyper.speciesIdentification.get_refseq_mash()
    ROOT_DIR = os.path.abspath(os.path.join(TEST_ROOT, '..'))
    MASHSTATSMETAFILE=os.path.join(TEST_ROOT+"/mash_refseq_meta.txt")
    MASHINFILE = os.path.join(ROOT_DIR, 'ectyper/Data/refseq.genomes.k21s1000.msh')
    ASSEMBLYREFSEQMETAFILE = os.path.join(ROOT_DIR, 'ectyper/Data/assembly_summary_refseq.txt')


    cmd = ["mash info -t " +  MASHINFILE   + " > " + MASHSTATSMETAFILE]
    print("File written to {}".format(MASHSTATSMETAFILE))
    subprocess.run(cmd, shell=True)
    mashsketchdatadf = pandas.read_csv(MASHSTATSMETAFILE,sep="\t")
    mashaccessions=[re.findall("(GCF_\d+)\.+",item)[0] for item in mashsketchdatadf.iloc[:,2].values.tolist()]
    LOG.info("Extracted {} MASH RefSeq accessions".format(len(mashaccessions)))


    genomeassemblystatrefseqsdf = pandas.read_csv(ASSEMBLYREFSEQMETAFILE, sep="\t", skiprows=1)
    #print(genomeassemblystatrefseqsdf.iloc[:, 0].values.tolist()[0:10])
    metaaccessionsrefseq = [re.findall("(GCF_\d+)\.+",item)[0]  for item  in genomeassemblystatrefseqsdf.iloc[:, 0].values.tolist()]
    #print(metaaccessionsrefseq[0:10])
    metaaccessionsrefseqdict = dict.fromkeys(metaaccessionsrefseq, True)


    LOG.info("Extracted {} assembly meta data accessions".format(len(metaaccessionsrefseq)))

    notfoundaccessions = []
    for accession in mashaccessions:
        if not metaaccessionsrefseqdict.get(accession):
            notfoundaccessions.append(accession)

    LOG.info("Not found accessions in meta data {}".format(len(notfoundaccessions)))
    print(notfoundaccessions)

    #genomeassemblystatsdf.loc[1:2, ["# assembly_accession",'bioproject', 'biosample',"wgs_master",
    # "refseq_category","taxid","species_taxid", "organism_name","ftp_path"]]



    # assembly_accession    bioproject      biosample       wgs_master      refseq_category taxid   species_taxid   organism_name   infraspecific_name      isolate version_status  assembly_level  release_type    genome_rep      seq_rel_date    asm_name        submitter       gbrs_paired_asm paired_asm_comp ftp_path        excluded_from_refseq    relation_to_type_material