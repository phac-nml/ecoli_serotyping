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


def test_getGenomeName():

    assert getGenomeName('AAJT0203.1|gi|190904618|gb|AAJT02000002.1|emb|AAJT0201|dbj|AAJT0202|lcl|AAJT0204|ref|AA_JT0205| Escherichia coli B7A gcontig_1112495745514, whole genome shotgun sequence', 'AAJT02') == 'AAJT0204'

    assert getGenomeName('AAJT0203.1|gi|190904618|gb|AAJT02000002.1|emb|AAJT0201|dbj|AAJT0202|ref|AA_JT0205| Escherichia coli B7A gcontig_1112495745514, whole genome shotgun sequence', 'AAJT02') == 'AAJT0203.1'

    assert getGenomeName('gi|190904618|gb|AAJT02000002.1|emb|AAJT0201|dbj|AAJT0202|ref|AA_JT0205| Escherichia coli B7A gcontig_1112495745514, whole genome shotgun sequence', 'AAJT02') == 'AAJT0200'

    assert getGenomeName('gi|190904618|emb|AAJT0201|dbj|AAJT0202|ref|AA_JT0205 Escherichia coli B7A gcontig_1112495745514, whole genome shotgun sequence', 'AAJT02') == 'AAJT0201'

    assert getGenomeName('gi|190904618|dbj|AAJT0202|ref|AA_JT0205| Escherichia coli B7A gcontig_1112495745514, whole genome shotgun sequence', 'AAJT02') == 'AAJT0202'

    assert getGenomeName('gi|190904618|ref|AA_JT0205| Escherichia coli B7A gcontig_1112495745514, whole genome shotgun sequence', 'AAJT02') == 'AA_JT0205'

    assert getGenomeName('gi|190904618| Escherichia coli B7A gcontig_1112495745514, whole genome shotgun sequence', 'AAJT02' ) == '19090461'

    assert getGenomeName('Escherichia coli B7A gcontig_1112495745514, whole genome shotgun sequence', 'AAJT02') == 'AAJT02'


def test_checkFiles():

    assert checkFiles(getListGenomes(REL_DIR)) == expectedList1

    with raises(SystemExit):
        checkFiles(getListGenomes(SCRIPT_DIRECTORY + '../Data/Testing_Data/Invalid Files/'))

    assert checkFiles(getListGenomes(SCRIPT_DIRECTORY + '../Data/Testing_Data')) == expectedList1



def test_initializeDB():

    os.chdir('/home/calarose/Documents/')

    assert initializeDB() == 0

    os.chdir('/home/calarose/ectyper/Data/')

    assert initializeDB() == 0

    os.chdir('/home/calarose/ectyper/')

    assert initializeDB() == 0

# def test_runBlastQuery():
#
# def test_parseResults():
#
# def test_findPerfectMatches():






