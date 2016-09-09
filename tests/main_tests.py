from ectyper.src.ectyper import *
import os

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
REL_DIR = '../Data/Reference Genomes/'

def test_getListGenomes():

    sortedExpectedList = sorted([
          SCRIPT_DIRECTORY + REL_DIR + 'NC_011749.1.fasta'
        , SCRIPT_DIRECTORY + REL_DIR + 'NC_011739.1.fasta'
        , SCRIPT_DIRECTORY + REL_DIR + 'NC_011751.1.fasta'
        , SCRIPT_DIRECTORY + REL_DIR + 'AAJT02.1.fsa_nt'
        , SCRIPT_DIRECTORY + REL_DIR + 'NC_002695.1.fasta'
        , SCRIPT_DIRECTORY + REL_DIR + 'NC_002128.1.fasta'
        , SCRIPT_DIRECTORY + REL_DIR + 'NC_011750.1.fasta'
        , SCRIPT_DIRECTORY + REL_DIR + 'NC_002127.1.fasta'])

    assert getListGenomes(SCRIPT_DIRECTORY + REL_DIR) == sortedExpectedList

    assert getListGenomes(SCRIPT_DIRECTORY + REL_DIR + 'NC_011750.1.fasta') == [os.path.abspath(SCRIPT_DIRECTORY + REL_DIR + 'NC_011750.1.fasta')]

    assert getListGenomes(SCRIPT_DIRECTORY +'../Data/') == sortedExpectedList

    assert getListGenomes(SCRIPT_DIRECTORY +'../Data/Empty Folder/') == []

# def test_checkFiles():
#     assert checkList() ==
#     assert checklist() ==
#     assert checklist() ==
#     assert checklist() ==









