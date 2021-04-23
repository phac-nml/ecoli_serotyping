import sys
import os
import logging
import re
from ectyper import ectyper
import tempfile
LOG=logging.getLogger("TEST")
TEST_ROOT = os.path.dirname(__file__)

def set_input(input,
              percent_iden=None,
              output=tempfile.mkdtemp(),
              customsketchfile=None,
              cores=1,
              print_sequence=False,
              verify=False,
              debug=False):
    """
    Create the sys.argv[] without need for commandline input.
    :param input: Input file given by testing function
    :param percent_iden: Percent identity for comparison
    :param output: Location of output
    :return: None
    """
    args = ['-i', input,
            '-c', str(cores)
            ]
    if verify:
        args+=['--verify']
    if percent_iden:
        args += ['-d', str(percent_iden)]
    if output:
        args += ['-o', output]
    if print_sequence:
        args += ['--sequence']
    if customsketchfile:
        args += ['-r', customsketchfile]
    if debug:
        args+=['--debug']
    sys.argv[1:] = args


def test_Otyping(caplog):
    """
    Giving E.coli fasta genomes with truncated wzx and wzy genes with reference coverage <50 predict O and H antigens
    :return: None
    """
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT, 'Data/Escherichia_O26H11.fasta')#+","+os.path.join(TEST_ROOT, 'Data/Escherichia.fna')
    tmpdir=tempfile.mkdtemp()
    set_input(input=file,cores=4,print_sequence=True, verify=True, output=tmpdir, debug=False)

    ectyper.run_program()


    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
         secondrow = outfp.readlines()[1].split("\t")
         Otype = secondrow[2]
         Htype = secondrow[3]

    assert Otype == "-", "Expected no call but reported O-type:" + Otype
    assert Htype == "H11", "Expected H11 but reported H-type:" + Htype


def test_closeOalles_O42_O28(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/EscherichiaO28H5.fasta')  # +","+os.path.join(TEST_ROOT, 'Data/Escherichia.fna')

    tmpdir = tempfile.mkdtemp()
    set_input(input=file, cores=4, print_sequence=True, verify=True, debug=False, output=tmpdir)
    ectyper.run_program()
    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
         secondrow = outfp.readlines()[1]
    print(secondrow)
    assert re.match(r".+Escherichia coli.+O42\/O28\tH25\tO42\/O28:H25", secondrow)


def test_Shigella_typing(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/DRR015915_Shigella_boydii.fasta')  # +","+os.path.join(TEST_ROOT, 'Data/Escherichia.fna')
    tmpdir = tempfile.mkdtemp()
    set_input(input=file, cores=4, print_sequence=True, debug=True, verify=True, output=tmpdir)
    ectyper.run_program()

    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
         secondrow = outfp.readlines()[1].split("\t")
         species = secondrow[1]
    assert species == "Shigella boydii"

def test_download_refseq_mash(caplog):
    caplog.set_level(logging.DEBUG)
    response = ectyper.speciesIdentification.get_refseq_mash_and_assembly_summary()
    assert response == True,"Something went wrong with the RefSeq sketch download. Check internet connection ..."


def test_mixofspecies(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/Campylobacter.fasta') +","+os.path.join(TEST_ROOT, 'Data/Salmonella.fasta')+","\
                         + os.path.join(TEST_ROOT, 'Data/Escherichia.fastq')
    tmpdir = tempfile.mkdtemp()
    set_input(input=file, cores=4, print_sequence=True, verify=True, output=tmpdir)

    ectyper.run_program()

    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
         rows = outfp.readlines()
    rows=rows[1:] #remove header line

    serovars=[]; genomenames=[]; QCflag=[]; confidence=[]
    for row in rows:
        rowlist = row.split("\t")
        print(rowlist)
        serovars.append(rowlist[4])
        genomenames.append(rowlist[1])
        QCflag.append(rowlist[5])
        confidence.append(rowlist[6])

    assert serovars == ['-:-', 'O22:H8', '-:-']
    expectedspecies_list = ["Campylobacter jejuni","Escherichia coli","Salmonella enterica"]
    for i in range(0,3):
        assert bool(re.match(expectedspecies_list[i], genomenames[i])) == True
    assert QCflag == ["WARNING (WRONG SPECIES)","PASS (REPORTABLE)","WARNING (WRONG SPECIES)"]


def test_Ealbertii_1(caplog): #error
    LOG.info("Starting 1 of 3 test on EnteroBase on sample ESC_HA8355AA_AS: Escherichia albertii O65:H5")
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

def test_Ealbertii_2(): #error
    LOG.info("Starting 2 of 3 test on EnteroBase on sample on ESC_HA8509AA_AS: Escherichia albertii O5:H5")

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
    LOG.info("Starting 3 of 3 test Escherichia albertii O49:NM") #can not type O49 due to poor sequence quality of uncertainty of wet-lab O49 typing
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


def test_Ecoli_O17H18(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/EscherichiaO17H18.fasta')
    tmpdir = tempfile.mkdtemp()
    set_input(input=file, cores=4, print_sequence=False, verify=True, debug=True, output=tmpdir)

    ectyper.run_program()

    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
         rows = outfp.readlines()
    secondrow=rows[1:][0] #check only second row
    assert "Escherichia coli\tO77/O17/O44/O106\tH18\tO77/O17/O44/O106:H18\tWARNING MIXED O-TYPE" in secondrow


