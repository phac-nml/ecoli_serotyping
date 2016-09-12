from pytest import raises

from ectyper.src.ectyper import *
import os

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
REL_DIR = SCRIPT_DIRECTORY + "../Data/Testing_Data/Reference Genomes"

expectedList1 = sorted([
    REL_DIR + '/NC_011749.1.fasta'
    , REL_DIR + '/NC_011739.1.fasta'
    , REL_DIR + '/NC_011751.1.fasta'
    , REL_DIR + '/AAJT02.1.fsa_nt'
    , REL_DIR + '/NC_002695.1.fasta'
    , REL_DIR + '/NC_002128.1.fasta'
    , REL_DIR + '/NC_011750.1.fasta'
    , REL_DIR + '/NC_002127.1.fasta'])

expectedList2 = sorted([
    SCRIPT_DIRECTORY + '../Data/Testing_Data/Invalid Files/Invalid-File-1.txt',
    SCRIPT_DIRECTORY + '../Data/Testing_Data/Invalid Files/Rabbit-in-Hat.jpg'
])


def test_getListGenomes():

    assert getListGenomes(REL_DIR) == expectedList1

    assert getListGenomes(REL_DIR + '/NC_011750.1.fasta') == [os.path.abspath(REL_DIR + '/NC_011750.1.fasta')]

    expectedList2.extend(expectedList1)

    assert getListGenomes(SCRIPT_DIRECTORY +'../Data/Testing_Data/') == sorted(expectedList2)

    assert getListGenomes(SCRIPT_DIRECTORY +'../Data/Testing_Data/Empty Folder/') == []


def test_checkFiles():

    assert checkFiles(getListGenomes(REL_DIR)) == expectedList1

    with raises(SystemExit):
        checkFiles(getListGenomes(SCRIPT_DIRECTORY + '../Data/Testing_Data/Invalid Files/'))

    assert checkFiles(getListGenomes(SCRIPT_DIRECTORY + '../Data/Testing_Data')) == expectedList1









