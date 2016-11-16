
import hashlib
import os.path
import pytest
import json

from src.Serotyper.ecvalidatingfiles import *

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
REL_DIR = SCRIPT_DIRECTORY + "../Data/Testing_Data/Reference_Genomes/"

expectedList1 = sorted([
    REL_DIR + 'NC_011749.1.fasta'
    , REL_DIR + 'NC_011739.1.fasta'
    , REL_DIR + 'NC_011751.1.fasta'
    , REL_DIR + 'AAJT02.1.fsa_nt'
    , REL_DIR + 'NC_002695.1.fasta'
    , REL_DIR + 'NC_002128.1.fasta'
    , REL_DIR + 'NC_011750.1.fasta'
    , REL_DIR + 'NC_002127.1.fasta'])

expectedList2 = sorted([
    SCRIPT_DIRECTORY + '../Data/Testing_Data/Invalid Files/Invalid-File-1.txt',
    SCRIPT_DIRECTORY + '../Data/Testing_Data/Invalid Files/Rabbit-in-Hat.jpg'
])

REL_DIR2 = SCRIPT_DIRECTORY + '../temp/xml/'

expectedList3 = sorted([os.path.abspath(REL_DIR2 + 'AAJT02.1.xml'),
                        os.path.abspath(REL_DIR2 + 'NC_002127.1.xml'),
                        os.path.abspath(REL_DIR2 + 'NC_002128.1.xml'),
                        os.path.abspath(REL_DIR2 + 'NC_002695.1.xml'),
                        os.path.abspath(REL_DIR2 + 'NC_011739.1.xml'),
                        os.path.abspath(REL_DIR2 + 'NC_011749.1.xml'),
                        os.path.abspath(REL_DIR2 + 'NC_011750.1.xml'),
                        os.path.abspath(REL_DIR2 + 'NC_011751.1.xml')])

def test_getFilesList():

    if not os.path.isdir(SCRIPT_DIRECTORY +'../Data/Testing_Data/Empty_Folder/'):
        os.mkdir(SCRIPT_DIRECTORY +'../Data/Testing_Data/Empty_Folder/')

    assert getFilesList(REL_DIR) == expectedList1

    assert getFilesList(REL_DIR + 'NC_011750.1.fasta') == [os.path.abspath(REL_DIR + '/NC_011750.1.fasta')]

    expectedList2.extend(expectedList1)

    assert getFilesList(SCRIPT_DIRECTORY +'../Data/Testing_Data/') == sorted(expectedList2)

    assert getFilesList(SCRIPT_DIRECTORY +'../Data/Testing_Data/Empty_Folder/') == []


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

    assert checkFiles(getFilesList(REL_DIR)) == expectedList1

    assert checkFiles(getFilesList(SCRIPT_DIRECTORY + '../Data/Testing_Data/Invalid Files/')) == 'Error'

    assert checkFiles(getFilesList(SCRIPT_DIRECTORY + '../Data/Testing_Data')) == expectedList1



def test_initializeDB():

    assert subprocess.call(["/usr/bin/makeblastdb", "-in", SCRIPT_DIRECTORY + "../Data/EcOH.fasta ", "-dbtype", "nucl", "-title", "ECTyperDB"]) == 1

    assert subprocess.call(["/usr/bin/makeblastdb", "-in", SCRIPT_DIRECTORY + "../Data/EcOH.fasta ", "-dbtype", "-title", "ECTyperDB", "-out", REL_DIR + "ECTyperDB"]) == 1

    assert subprocess.call(["/usr/bin/makeblastdb", "-in", SCRIPT_DIRECTORY , "-dbtype", "nucl", "-title", "ECTyperDB", "-out", REL_DIR + "ECTyperDB"]) == 1


    # os.chdir('/home/calarose/Documents/')
    # assert initializeDB() == 0
    #
    # os.chdir('/home/calarose/ectyper/Data/')
    # assert initializeDB() == 0
    #
    # os.chdir('/home/calarose/ectyper/')
    # assert initializeDB() == 0


def test_runBlastQuery():

   if not os.path.isdir(SCRIPT_DIRECTORY+ "../temp/"):
        os.mkdir(SCRIPT_DIRECTORY+ "../temp/")
   if not os.path.isdir(SCRIPT_DIRECTORY + "../temp/xml/"):
       os.mkdir(SCRIPT_DIRECTORY + "../temp/xml/")
   if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/Uploads/'):
       os.mkdir(SCRIPT_DIRECTORY + '../temp/Uploads/')

   initializeDB()
   assert runBlastQuery(expectedList1, 'ECTyperDB') == os.path.abspath(SCRIPT_DIRECTORY + '../temp/xml/combined_genomes.xml')
   assert runBlastQuery([REL_DIR + '/AAJT02.1.fsa_nt'], 'ECTyperDB') == os.path.abspath(SCRIPT_DIRECTORY + '../temp/xml/AAJT02.1.xml')


def test_parseResults():

    with open(SCRIPT_DIRECTORY + '../Data/Test_dictionaries/ecvalidatingfiles_dict.json') as f:
        expectedDict = json.load(f)

    testDict = parseResults(runBlastQuery([REL_DIR + 'NC_011749.1.fasta',REL_DIR + 'AAJT02.1.fsa_nt'], 'ECTyperDB'))

    for test_genome, test_alignment in testDict.iteritems():
        if test_genome in expectedDict:
            if isinstance(test_alignment, dict):
                for test_title, test_info in test_alignment.iteritems():
                    test_length = test_info.keys()[0]
                    test_hsp = test_info[test_length]
                    hash_hsp = hashlib.md5()
                    hash_hsp.update(str(test_hsp))
                    assert expectedDict[test_genome][test_title] == str(hash_hsp.hexdigest())
                    # if expectedDict[test_genome][test_title] == str(hash_hsp.hexdigest()):
                    #     pytest.fail("The dictionaries aren't the same. Test failed. \nCAUSE \nAlignment title: " + str(test_title)  + '\nExpected: ' + str(hash_hsp.hexdigest()) +  '\nFound: ' +
                    #                 str(expectedDict[test_genome][test_title]) )
        else:
            pytest.fail("The dictionaries aren't the same. Test failed. \nCAUSE \n" + str(test_genome) + " isn't a key")







