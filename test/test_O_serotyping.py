import sys
import os
import logging
import re
from ectyper import ectyper, definitions
import tempfile
LOG=logging.getLogger("TEST")
LOG.setLevel(logging.INFO)
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
         lines = outfp.readlines()
         header = lines[0].split("\t")
         secondrow = lines[1].split("\t")
         assert 'O-type' in header, f'O-type column not found in the output.tsv ({header})'
         assert 'H-type' in header, f'H-type column not found in the output.tsv ({header})'
         col_number = [idx for idx, i in enumerate(header) if i == "O-type"][0]
         Otype = secondrow[col_number]
         col_number = [idx for idx, i in enumerate(header) if i == "H-type"][0]
         Htype = secondrow[col_number]

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
    assert re.match(r".+Escherichia coli.+O28\/O42\tH25\tO28\/O42:H25", secondrow)


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
    header = rows[0]
    resultsrows = rows[1:] #remove header line

    serotypes=[]; species=[]; QCflag=[]
         
    for row in resultsrows:
        rowlist = row.split("\t")
        serotypes.append(rowlist[ [idx for idx, i in enumerate(header.split('\t')) if i == "Serotype"][0] ])
        species.append(rowlist[  [idx for idx, i in enumerate(header.split('\t')) if i == "Species"][0]   ])
        QCflag.append(rowlist[ [idx for idx, i in enumerate(header.split('\t')) if i == "QC"][0]  ])

    assert serotypes == ['-:-', 'O22:H8', '-:-']
    expectedspecies_list = ["Campylobacter_D jejuni","Escherichia coli","Salmonella enterica"]
    for i in range(0,3):
        assert bool(re.match(expectedspecies_list[i], species[i])) == True
    assert QCflag == ["WARNING (WRONG SPECIES)","PASS (REPORTABLE)","WARNING (WRONG SPECIES)"]

def test_Ecoli_O17H18(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/EscherichiaO17H18.fasta')
    tmpdir = tempfile.mkdtemp()
    set_input(input=file, cores=4, print_sequence=False, verify=True, debug=True, output=tmpdir)

    ectyper.run_program()

    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
         rows = outfp.readlines()
    columnslist = rows[0].split('\t')
    secondrow=rows[1] #check only second row
    secondrow_values=secondrow.split('\t')
    assert len(secondrow_values) == len(columnslist), f"Unequal number of columns {len(columnslist)} and values {len(secondrow_values)}"
    assert "Escherichia coli" in secondrow.split('\t')
    assert "O17/O44/O77/O106\tH18\tO17/O44/O77/O106:H18\tWARNING MIXED O-TYPE" in secondrow

def test_Ecoli_O178H19_singleOant_highsimilarity(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/EC20121408_O178H19_singleO178.fa.gz') #EC20180127_O153O178H19.fa
    tmpdir = tempfile.mkdtemp()
    set_input(input=file, cores=4, print_sequence=False, verify=True, debug=True, output=tmpdir)

    ectyper.run_program()

    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
         rows = outfp.readlines()
    
    secondrow=rows[1] #check only second row
    secondrow_values=secondrow.split('\t')
    assert "Escherichia coli" in secondrow_values
    assert "O178:H19" in secondrow_values
    
def test_Ecoli_O153O178H19_mixedOant_highsimilarity(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/EC20180127_O153O178H19_mixedO.fa.gz') 
    tmpdir = tempfile.mkdtemp()
    set_input(input=file, cores=4, print_sequence=False, verify=True, debug=True, output=tmpdir)

    ectyper.run_program()

    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
         rows = outfp.readlines()
    
    secondrow=rows[1] #check only second row
    secondrow_values=secondrow.split('\t')
    assert "Escherichia coli" in secondrow_values
    assert "O153/O178:H19" in secondrow_values



def test_download_refseq_mash(caplog, tmpdir):
    caplog.set_level(logging.DEBUG)
    response = ectyper.speciesIdentification.get_species_mash(definitions.SPECIES_ID_SKETCH)
    assert response == True,"Something went wrong with the Species ID database download. Check internet connection ..."    


