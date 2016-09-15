from pytest import raises
import hashlib
import pytest
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

REL_DIR2 = SCRIPT_DIRECTORY + '../temp/'
expectedList3 = sorted([os.path.abspath(REL_DIR2 + 'AAJT02.1.xml'),
                        os.path.abspath(REL_DIR2 + 'NC_002127.1.xml'),
                        os.path.abspath(REL_DIR2 + 'NC_002128.1.xml'),
                        os.path.abspath(REL_DIR2 + 'NC_002695.1.xml'),
                        os.path.abspath(REL_DIR2 + 'NC_011739.1.xml'),
                        os.path.abspath(REL_DIR2 + 'NC_011749.1.xml'),
                        os.path.abspath(REL_DIR2 + 'NC_011750.1.xml'),
                        os.path.abspath(REL_DIR2 + 'NC_011751.1.xml')])

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

    assert subprocess.call(["/usr/bin/makeblastdb", "-in", SCRIPT_DIRECTORY + "../Data/EcOH.fasta ", "-dbtype", "nucl", "-title", "ECTyperDB"]) == 1

    assert subprocess.call(["/usr/bin/makeblastdb", "-in", SCRIPT_DIRECTORY + "../Data/EcOH.fasta ", "-dbtype", "-title", "ECTyperDB", "-out", REL_DIR + "ECTyperDB"]) == 1

    assert subprocess.call(["/usr/bin/makeblastdb", "-in", SCRIPT_DIRECTORY , "-dbtype", "nucl", "-title", "ECTyperDB", "-out", REL_DIR + "ECTyperDB"]) == 1

    os.chdir('/home/calarose/Documents/')
    assert initializeDB() == 0

    os.chdir('/home/calarose/ectyper/Data/')
    assert initializeDB() == 0

    os.chdir('/home/calarose/ectyper/')
    assert initializeDB() == 0


def test_runBlastQuery():

   assert runBlastQuery(expectedList1) == expectedList3


def test_parseResults():

    expectedDict = {
        'gnl|BL_ORD_ID|175 8__wzx__wzx-O148__176 DQ167407.1;O antigen flippase;O148' : 'aee07ddae92b95074a93148a06624a29', 'gnl|BL_ORD_ID|377 9__wzy__wzy-O148__378 DQ167407.1;O antigen polyermase;O148' : 'c10c9586ec7845c902fbbd55f7df991f',
        'gnl|BL_ORD_ID|190 8__wzx__wzx-O159__191 EU294176.1;O antigen flippase;O159' : '8cab52763f6262bcecaa41bece2ae41f', 'gnl|BL_ORD_ID|165 8__wzx__wzx-O140__166 AB812060.1;O antigen flippase;O140' : '27d4fd2e1144f746437c3fc4bfcf7ec3',
        'gnl|BL_ORD_ID|28 1__fliC__fliC-H28__29 AY250010.1;flagellin;H28' : '6dfdab8081e7d834f1b3243985472b2d', 'gnl|BL_ORD_ID|31 1__fliC__fliC-H31__32 AY250013.1;flagellin;H31' : '0ef43266b0331ed986db7371dabdfc98',
        'gnl|BL_ORD_ID|10 1__fliC__fliC-H14__11 AY249998.1;flagellin;H14' : 'd4fb30f63599ca22f5ce3c6b9f8028d5', 'gnl|BL_ORD_ID|7 1__fliC__fliC-H12__8 AY337471.1;flagellin;H12' : 'b8c9c7ae76ae13eea5287cdd654d3907',
        'gnl|BL_ORD_ID|8 1__fliC__fliC-H12__9 AY337474.1;flagellin;H12' : 'ac09ce2077f2225059bc306c8c83318d', 'gnl|BL_ORD_ID|62 1__fliC__fliC-H6__63 AY249991.1;flagellin;H6' : 'a52719aff71d6a30d58d0e6a2023a995',
        'gnl|BL_ORD_ID|35 1__fliC__fliC-H34__36 AF34585;flagellin;H34' : 'eadeb19641cb38052520d6c92037790b', 'gnl|BL_ORD_ID|47 1__fliC__fliC-H41__48 AY250020.1;flagellin;H41' : '1f3cef59dafb2c38024d508ac27439fd',
        'gnl|BL_ORD_ID|34 1__fliC__fliC-H34__35 AY250016.1;flagellin;H34' : 'eadeb19641cb38052520d6c92037790b','gnl|BL_ORD_ID|9 1__fliC__fliC-H12__10 AY249997.1;flagellin;H12' : '779f01882e3901826062c3d55b8d319a',
        'gnl|BL_ORD_ID|0 1__fliC__fliC-H1__1 AB028471.1;flagellin;H1' : '6143271ae1219405eab366766d4b7593', 'gnl|BL_ORD_ID|50 1__fliC__fliC-H45__51 AY250023.1;flagellin;H45' : 'ad67533120333388945a7712ea6a53b2',
        'gnl|BL_ORD_ID|1 1__fliC__fliC-H1__2 L07387.1;flagellin;H1' : '61dc9c31e66819e245cb7a9ca2e70f1d', 'gnl|BL_ORD_ID|55 1__fliC__fliC-H49__56 AY250026.1;flagellin;H49' : '8a06559dd8a39a175ed0aaaf58e83c6f',
        'gnl|BL_ORD_ID|54 1__fliC__fliC-H49__55 AB028480.1;flagellin;H49' : '8a06559dd8a39a175ed0aaaf58e83c6f', 'gnl|BL_ORD_ID|52 1__fliC__fliC-H46__53 AY250024.1;flagellin;H46' : '780b60fa959c1550fd8bdb0f39b5f579',
        'gnl|BL_ORD_ID|66 1__fliC__fliC-H7__67 AF228494.1;flagellin;H7' : '3355cec6378330c6530aedc2700af69e', 'gnl|BL_ORD_ID|64 1__fliC__fliC-H7__65 AF228491.1;flagellin;H7' : '3355cec6378330c6530aedc2700af69e',
        'gnl|BL_ORD_ID|75 1__fliC__fliC-H9__76 AY249994.1;flagellin;H9' : '1489feeac76fb141ac77e2163ecbf1a9', 'gnl|BL_ORD_ID|51 1__fliC__fliC-H46__52 AB028478.1;flagellin;H46' : 'bbcc42bbde8a226400ae3494b900880a',
        'gnl|BL_ORD_ID|72 1__fliC__fliC-H7__73 AY249992.1;flagellin;H7' : '86d2704ec42a318a7b26528261465ff5','gnl|BL_ORD_ID|65 1__fliC__fliC-H7__66 AF228492.1;flagellin;H7' : '86d2704ec42a318a7b26528261465ff5',
        'gnl|BL_ORD_ID|11 1__fliC__fliC-H15__12 AY249999.1;flagellin;H15' : 'cb457fe7fb890888e83ca1432af75528', 'gnl|BL_ORD_ID|68 1__fliC__fliC-H7__69 AF228496.1;flagellin;H7' : 'c8fd3b49aaadbef6c114577eb35ba872',
        'gnl|BL_ORD_ID|60 1__fliC__fliC-H52__61 AY250028.1;flagellin;H52' : 'a04d1a965df7932b0b42fa10764ec3d0', 'gnl|BL_ORD_ID|71 1__fliC__fliC-H7__72 AB334575.1;flagellin;H7' : '3fb9b7893b507725c49c4be19002b831',
        'gnl|BL_ORD_ID|70 1__fliC__fliC-H7__71 AB334574.1;flagellin;H7' : '3fb9b7893b507725c49c4be19002b831', 'gnl|BL_ORD_ID|69 1__fliC__fliC-H7__70 AY337468.1;flagellin;H7' : '3fb9b7893b507725c49c4be19002b831',
        'gnl|BL_ORD_ID|67 1__fliC__fliC-H7__68 AF228495.1;flagellin;H7' : '3fb9b7893b507725c49c4be19002b831', 'gnl|BL_ORD_ID|63 1__fliC__fliC-H7__64 AF228487.1;flagellin;H7' : '3fb9b7893b507725c49c4be19002b831',
        'gnl|BL_ORD_ID|26 1__fliC__fliC-H26__27 AY250008.1;flagellin;H26' : '11bec9ac7863945e285044b85fa328bc', 'gnl|BL_ORD_ID|21 1__fliC__fliC-H20__22 AY250003.1;flagellin;H20' : 'c236f038a83412f0a51e6776f6ca2b9b',
        'gnl|BL_ORD_ID|16 1__fliC__fliC-H18__17 AY250001.1;flagellin;H18' : 'a16b6e1c69abbb85dcc3c43d5cc28552', 'gnl|BL_ORD_ID|18 1__fliC__fliC-H19__19 AY250002.1;flagellin;H19' : '32eee96c2e861bfd9c1876647b3d5d3e',
        'gnl|BL_ORD_ID|17 1__fliC__fliC-H19__18 AY337479.1;flagellin;H19' : '32eee96c2e861bfd9c1876647b3d5d3e', 'gnl|BL_ORD_ID|36 1__fliC__fliC-H37__37 AY250017.1;flagellin;H37' : 'bd80deb2fdbfd9d46ea82408be015218',
        'gnl|BL_ORD_ID|23 1__fliC__fliC-H23__24 AY250005.1;flagellin;H23' : 'eff1f961556b4a88b8072371521275c8', 'gnl|BL_ORD_ID|30 1__fliC__fliC-H30__31 AY250011.1;flagellin;H30' : 'e02e1735623030f5952d26e0f6d3b90d',
        'gnl|BL_ORD_ID|32 1__fliC__fliC-H32__33 AY250014.1;flagellin;H32' : 'ed23dd977c2e0f3cbbacd066e1a3e39c', 'gnl|BL_ORD_ID|53 1__fliC__fliC-H48__54 AY250025.1;flagellin;H48' : 'aa7a0638df03c38e5db1dcb09d7fa3a3',
        'gnl|BL_ORD_ID|49 1__fliC__fliC-H43__50 AY250022.1;flagellin;H43' : '7c0cab8d80f814a75e241cd856aa22d7', 'gnl|BL_ORD_ID|59 1__fliC__fliC-H51__60 AY250027.1;flagellin;H51' : '1dddbe105cd8f01595d077fb3e96ab1d',
        'gnl|BL_ORD_ID|58 1__fliC__fliC-H51__59 AB028481.1;flagellin;H51' : '8568a7d4f80bbe9671256f44107a713c', 'gnl|BL_ORD_ID|81 3__fllA__fllA-H44__82 #VALUE!;flagellin;H44' : '50034b78fa6139f4049455b3bf35d759',
        'gnl|BL_ORD_ID|3 1__fliC__fliC-H10__4 AY249995.1;flagellin;H10' : '46d509f1c57c5cb74e9edf1a41760238', 'gnl|BL_ORD_ID|2 1__fliC__fliC-H10__3 AY337482.1;flagellin;H10' : '46d509f1c57c5cb74e9edf1a41760238',
        'gnl|BL_ORD_ID|82 3__fllA__fllA-H55__83 AB269771;flagellin;H55' : '915d2ce61cdfab6a05054b48c479f306', 'gnl|BL_ORD_ID|57 1__fliC__fliC-H5__58 AY249990.1;flagellin;H5' : '2097d002ad1b83c687e03d7004d846d6',
        'gnl|BL_ORD_ID|61 1__fliC__fliC-H56__62 AY250029.1;flagellin;H56' : '56cab81dfe4a766531738f8747bb0caa', 'gnl|BL_ORD_ID|37 1__fliC__fliC-H38__38 AY337466.1;flagellin;H38' : 'fa132af24eb3d853f0b92d4da688a981',
        'gnl|BL_ORD_ID|39 1__fliC__fliC-H38__40 AY250018.1;flagellin;H38' : '94ab1c21417d8878ae8e3b0de44c88f5', 'gnl|BL_ORD_ID|38 1__fliC__fliC-H38__39 AY337473.1;flagellin;H38' : '94ab1c21417d8878ae8e3b0de44c88f5',
        'gnl|BL_ORD_ID|44 1__fliC__fliC-H4__45 AY249989.1;flagellin;H4' : 'cd3369f667072a2d8d54f8562ed29cc3', 'gnl|BL_ORD_ID|42 1__fliC__fliC-H4__43 AJ605764.1;flagellin;H4' : 'cd3369f667072a2d8d54f8562ed29cc3',
        'gnl|BL_ORD_ID|41 1__fliC__fliC-H4__42 AJ536600.1;flagellin;H4' : 'cd3369f667072a2d8d54f8562ed29cc3', 'gnl|BL_ORD_ID|56 1__fliC__fliC-H5__57 AY337469.1;flagellin;H5' : '9e28653f90c9e14104892dd3e94c5456',
        'gnl|BL_ORD_ID|43 1__fliC__fliC-H4__44 AJ605765.1;flagellin;H4' : 'f11f18caebb5c00b2a743df77b9c30a6', 'gnl|BL_ORD_ID|33 1__fliC__fliC-H33__34 AY250015.1;flagellin;H33' : 'c2d5dbb952d454e42878c16a49569a4d',
        'gnl|BL_ORD_ID|48 1__fliC__fliC-H42__49 AY250021.1;flagellin;H42' : 'd877438aa5b2ce1b83a29ee661c33849', 'gnl|BL_ORD_ID|29 1__fliC__fliC-H29__30 AY250012.1;flagellin;H29' : 'd877438aa5b2ce1b83a29ee661c33849',
        'gnl|BL_ORD_ID|84 5__flnA__flnA-H17__85 CP002291;flagellin;H17' : '0ba4eb36a2fcd2d9333cae1936431eec', 'gnl|BL_ORD_ID|40 1__fliC__fliC-H39__41 AY250019.1;flagellin;H39' : 'a45424bf04fc3aa95036e92b25b7a276',
        'gnl|BL_ORD_ID|45 1__fliC__fliC-H40__46 AJ865464.1;flagellin;H40' : 'd1a780bc0db40b1a4c90c7cb8ee32001', 'gnl|BL_ORD_ID|24 1__fliC__fliC-H24__25 AY250006.1;flagellin;H24' : '5282b0acfb371dc93c2e41fee4a785f1',
        'gnl|BL_ORD_ID|73 1__fliC__fliC-H8__74 AJ865465.1;flagellin;H8' : '9754f3d51f8129ec865f5c6828ab209e', 'gnl|BL_ORD_ID|5 1__fliC__fliC-H11__6 AY337472.1;flagellin;H11' : '8caff7bb15fd57658b9edebdd18d6c2e',
        'gnl|BL_ORD_ID|4 1__fliC__fliC-H11__5 AY337465.1;flagellin;H11' : '8caff7bb15fd57658b9edebdd18d6c2e', 'gnl|BL_ORD_ID|74 1__fliC__fliC-H8__75 AJ884569.1;flagellin;H8' : '86b48ff8c3149dc361923c8b89a23563',
        'gnl|BL_ORD_ID|46 1__fliC__fliC-H40__47 AJ884568.1;flagellin;H40' : '655ca62dbbbdaac23f26bc54c1f51659', 'gnl|BL_ORD_ID|20 1__fliC__fliC-H2__21 HQ116826;flagellin;H2' : '8237e98abe46896f7bb6d8ddeb3275e1',
        'gnl|BL_ORD_ID|19 1__fliC__fliC-H2__20 AF543692;flagellin;H2' : 'c4966dee16a320fc55a3f9619b74dc12', 'gnl|BL_ORD_ID|27 1__fliC__fliC-H27__28 AY250009.1;flagellin;H27' : 'f3939ff5d0e1231ae3d89ed15e3e9d25',
        'gnl|BL_ORD_ID|22 1__fliC__fliC-H21__23 AY250004.1;flagellin;H21' : '8b4b3b29ab37a007ca04a57923520032', 'gnl|BL_ORD_ID|6 1__fliC__fliC-H11__7 AY249996.1;flagellin;H11' : '8b4b3b29ab37a007ca04a57923520032',
        'gnl|BL_ORD_ID|14 1__fliC__fliC-H16__15 AY337477.1;flagellin;H16' : '09064bc24643dc9bb5ffa42624ec78d5', 'gnl|BL_ORD_ID|13 1__fliC__fliC-H16__14 AY337476.1;flagellin;H16' : '09064bc24643dc9bb5ffa42624ec78d5' ,
        'gnl|BL_ORD_ID|12 1__fliC__fliC-H16__13 AY337475.1;flagellin;H16' : '09064bc24643dc9bb5ffa42624ec78d5', 'gnl|BL_ORD_ID|15 1__fliC__fliC-H16__16 AY250000.1;flagellin;H16' : '51f3c3538e86bd56fd16761505da4630',
        'gnl|BL_ORD_ID|25 1__fliC__fliC-H25__26 AY250007.1;flagellin;H25' : '94f0f51796d70c5e0a890beb107786e0' }

    testDict = parseResults([REL_DIR2 + 'NC_011749.1.xml',REL_DIR2 + 'AAJT02.1.xml'])

    for test_title, test_hsp in testDict.iteritems():
        hash_hsp = hashlib.md5()
        hash_hsp.update(str(test_hsp))
        if expectedDict[test_title] != str(hash_hsp.hexdigest()):
            pytest.fail("The dictionaries aren't the same. Test failed. \nCAUSE \nAlignment title: " + str(test_title) + "\nHSP:\n" + str(test_hsp))

def test_findPerfectMatches():

    expected_identical = {'gnl|BL_ORD_ID|377 9__wzy__wzy-O148__378 DQ167407.1;O antigen polyermase;O148' : 'c10c9586ec7845c902fbbd55f7df991f',
                     'gnl|BL_ORD_ID|190 8__wzx__wzx-O159__191 EU294176.1;O antigen flippase;O159' : '8cab52763f6262bcecaa41bece2ae41f' }

    expected_prediction = {'gnl|BL_ORD_ID|175 8__wzx__wzx-O148__176 DQ167407.1;O antigen flippase;O148' : 'aee07ddae92b95074a93148a06624a29', 'gnl|BL_ORD_ID|165 8__wzx__wzx-O140__166 AB812060.1;O antigen flippase;O140' : '27d4fd2e1144f746437c3fc4bfcf7ec3',
                       'gnl|BL_ORD_ID|28 1__fliC__fliC-H28__29 AY250010.1;flagellin;H28' : '6dfdab8081e7d834f1b3243985472b2d', 'gnl|BL_ORD_ID|31 1__fliC__fliC-H31__32 AY250013.1;flagellin;H31' : '0ef43266b0331ed986db7371dabdfc98',
                       'gnl|BL_ORD_ID|10 1__fliC__fliC-H14__11 AY249998.1;flagellin;H14' : 'd4fb30f63599ca22f5ce3c6b9f8028d5', 'gnl|BL_ORD_ID|7 1__fliC__fliC-H12__8 AY337471.1;flagellin;H12' : 'b8c9c7ae76ae13eea5287cdd654d3907',
                       'gnl|BL_ORD_ID|8 1__fliC__fliC-H12__9 AY337474.1;flagellin;H12' : 'ac09ce2077f2225059bc306c8c83318d', 'gnl|BL_ORD_ID|62 1__fliC__fliC-H6__63 AY249991.1;flagellin;H6' : 'a52719aff71d6a30d58d0e6a2023a995',
                       'gnl|BL_ORD_ID|35 1__fliC__fliC-H34__36 AF34585;flagellin;H34' : 'eadeb19641cb38052520d6c92037790b', 'gnl|BL_ORD_ID|47 1__fliC__fliC-H41__48 AY250020.1;flagellin;H41' : '1f3cef59dafb2c38024d508ac27439fd',
                       'gnl|BL_ORD_ID|34 1__fliC__fliC-H34__35 AY250016.1;flagellin;H34' : 'eadeb19641cb38052520d6c92037790b','gnl|BL_ORD_ID|9 1__fliC__fliC-H12__10 AY249997.1;flagellin;H12' : '779f01882e3901826062c3d55b8d319a',
                       'gnl|BL_ORD_ID|0 1__fliC__fliC-H1__1 AB028471.1;flagellin;H1' : '6143271ae1219405eab366766d4b7593', 'gnl|BL_ORD_ID|50 1__fliC__fliC-H45__51 AY250023.1;flagellin;H45' : 'ad67533120333388945a7712ea6a53b2',
                       'gnl|BL_ORD_ID|1 1__fliC__fliC-H1__2 L07387.1;flagellin;H1' : '61dc9c31e66819e245cb7a9ca2e70f1d', 'gnl|BL_ORD_ID|55 1__fliC__fliC-H49__56 AY250026.1;flagellin;H49' : '8a06559dd8a39a175ed0aaaf58e83c6f',
                       'gnl|BL_ORD_ID|54 1__fliC__fliC-H49__55 AB028480.1;flagellin;H49' : '8a06559dd8a39a175ed0aaaf58e83c6f', 'gnl|BL_ORD_ID|52 1__fliC__fliC-H46__53 AY250024.1;flagellin;H46' : '780b60fa959c1550fd8bdb0f39b5f579',
                       'gnl|BL_ORD_ID|66 1__fliC__fliC-H7__67 AF228494.1;flagellin;H7' : '3355cec6378330c6530aedc2700af69e', 'gnl|BL_ORD_ID|64 1__fliC__fliC-H7__65 AF228491.1;flagellin;H7' : '3355cec6378330c6530aedc2700af69e',
                       'gnl|BL_ORD_ID|75 1__fliC__fliC-H9__76 AY249994.1;flagellin;H9' : '1489feeac76fb141ac77e2163ecbf1a9', 'gnl|BL_ORD_ID|51 1__fliC__fliC-H46__52 AB028478.1;flagellin;H46' : 'bbcc42bbde8a226400ae3494b900880a',
                       'gnl|BL_ORD_ID|72 1__fliC__fliC-H7__73 AY249992.1;flagellin;H7' : '86d2704ec42a318a7b26528261465ff5','gnl|BL_ORD_ID|65 1__fliC__fliC-H7__66 AF228492.1;flagellin;H7' : '86d2704ec42a318a7b26528261465ff5',
                       'gnl|BL_ORD_ID|11 1__fliC__fliC-H15__12 AY249999.1;flagellin;H15' : 'cb457fe7fb890888e83ca1432af75528', 'gnl|BL_ORD_ID|68 1__fliC__fliC-H7__69 AF228496.1;flagellin;H7' : 'c8fd3b49aaadbef6c114577eb35ba872',
                       'gnl|BL_ORD_ID|60 1__fliC__fliC-H52__61 AY250028.1;flagellin;H52' : 'a04d1a965df7932b0b42fa10764ec3d0', 'gnl|BL_ORD_ID|71 1__fliC__fliC-H7__72 AB334575.1;flagellin;H7' : '3fb9b7893b507725c49c4be19002b831',
                       'gnl|BL_ORD_ID|70 1__fliC__fliC-H7__71 AB334574.1;flagellin;H7' : '3fb9b7893b507725c49c4be19002b831', 'gnl|BL_ORD_ID|69 1__fliC__fliC-H7__70 AY337468.1;flagellin;H7' : '3fb9b7893b507725c49c4be19002b831',
                       'gnl|BL_ORD_ID|67 1__fliC__fliC-H7__68 AF228495.1;flagellin;H7' : '3fb9b7893b507725c49c4be19002b831', 'gnl|BL_ORD_ID|63 1__fliC__fliC-H7__64 AF228487.1;flagellin;H7' : '3fb9b7893b507725c49c4be19002b831',
                       'gnl|BL_ORD_ID|26 1__fliC__fliC-H26__27 AY250008.1;flagellin;H26' : '11bec9ac7863945e285044b85fa328bc', 'gnl|BL_ORD_ID|21 1__fliC__fliC-H20__22 AY250003.1;flagellin;H20' : 'c236f038a83412f0a51e6776f6ca2b9b',
                       'gnl|BL_ORD_ID|16 1__fliC__fliC-H18__17 AY250001.1;flagellin;H18' : 'a16b6e1c69abbb85dcc3c43d5cc28552', 'gnl|BL_ORD_ID|18 1__fliC__fliC-H19__19 AY250002.1;flagellin;H19' : '32eee96c2e861bfd9c1876647b3d5d3e',
                       'gnl|BL_ORD_ID|17 1__fliC__fliC-H19__18 AY337479.1;flagellin;H19' : '32eee96c2e861bfd9c1876647b3d5d3e', 'gnl|BL_ORD_ID|36 1__fliC__fliC-H37__37 AY250017.1;flagellin;H37' : 'bd80deb2fdbfd9d46ea82408be015218',
                       'gnl|BL_ORD_ID|23 1__fliC__fliC-H23__24 AY250005.1;flagellin;H23' : 'eff1f961556b4a88b8072371521275c8', 'gnl|BL_ORD_ID|30 1__fliC__fliC-H30__31 AY250011.1;flagellin;H30' : 'e02e1735623030f5952d26e0f6d3b90d',
                       'gnl|BL_ORD_ID|32 1__fliC__fliC-H32__33 AY250014.1;flagellin;H32' : 'ed23dd977c2e0f3cbbacd066e1a3e39c', 'gnl|BL_ORD_ID|53 1__fliC__fliC-H48__54 AY250025.1;flagellin;H48' : 'aa7a0638df03c38e5db1dcb09d7fa3a3',
                       'gnl|BL_ORD_ID|49 1__fliC__fliC-H43__50 AY250022.1;flagellin;H43' : '7c0cab8d80f814a75e241cd856aa22d7', 'gnl|BL_ORD_ID|59 1__fliC__fliC-H51__60 AY250027.1;flagellin;H51' : '1dddbe105cd8f01595d077fb3e96ab1d',
                       'gnl|BL_ORD_ID|58 1__fliC__fliC-H51__59 AB028481.1;flagellin;H51' : '8568a7d4f80bbe9671256f44107a713c', 'gnl|BL_ORD_ID|81 3__fllA__fllA-H44__82 #VALUE!;flagellin;H44' : '50034b78fa6139f4049455b3bf35d759',
                       'gnl|BL_ORD_ID|3 1__fliC__fliC-H10__4 AY249995.1;flagellin;H10' : '46d509f1c57c5cb74e9edf1a41760238', 'gnl|BL_ORD_ID|2 1__fliC__fliC-H10__3 AY337482.1;flagellin;H10' : '46d509f1c57c5cb74e9edf1a41760238',
                       'gnl|BL_ORD_ID|82 3__fllA__fllA-H55__83 AB269771;flagellin;H55' : '915d2ce61cdfab6a05054b48c479f306', 'gnl|BL_ORD_ID|57 1__fliC__fliC-H5__58 AY249990.1;flagellin;H5' : '2097d002ad1b83c687e03d7004d846d6',
                       'gnl|BL_ORD_ID|61 1__fliC__fliC-H56__62 AY250029.1;flagellin;H56' : '56cab81dfe4a766531738f8747bb0caa', 'gnl|BL_ORD_ID|37 1__fliC__fliC-H38__38 AY337466.1;flagellin;H38' : 'fa132af24eb3d853f0b92d4da688a981',
                       'gnl|BL_ORD_ID|39 1__fliC__fliC-H38__40 AY250018.1;flagellin;H38' : '94ab1c21417d8878ae8e3b0de44c88f5', 'gnl|BL_ORD_ID|38 1__fliC__fliC-H38__39 AY337473.1;flagellin;H38' : '94ab1c21417d8878ae8e3b0de44c88f5',
                       'gnl|BL_ORD_ID|44 1__fliC__fliC-H4__45 AY249989.1;flagellin;H4' : 'cd3369f667072a2d8d54f8562ed29cc3', 'gnl|BL_ORD_ID|42 1__fliC__fliC-H4__43 AJ605764.1;flagellin;H4' : 'cd3369f667072a2d8d54f8562ed29cc3',
                       'gnl|BL_ORD_ID|41 1__fliC__fliC-H4__42 AJ536600.1;flagellin;H4' : 'cd3369f667072a2d8d54f8562ed29cc3', 'gnl|BL_ORD_ID|56 1__fliC__fliC-H5__57 AY337469.1;flagellin;H5' : '9e28653f90c9e14104892dd3e94c5456',
                       'gnl|BL_ORD_ID|43 1__fliC__fliC-H4__44 AJ605765.1;flagellin;H4' : 'f11f18caebb5c00b2a743df77b9c30a6', 'gnl|BL_ORD_ID|33 1__fliC__fliC-H33__34 AY250015.1;flagellin;H33' : 'c2d5dbb952d454e42878c16a49569a4d',
                       'gnl|BL_ORD_ID|48 1__fliC__fliC-H42__49 AY250021.1;flagellin;H42' : 'd877438aa5b2ce1b83a29ee661c33849', 'gnl|BL_ORD_ID|29 1__fliC__fliC-H29__30 AY250012.1;flagellin;H29' : 'd877438aa5b2ce1b83a29ee661c33849',
                       'gnl|BL_ORD_ID|84 5__flnA__flnA-H17__85 CP002291;flagellin;H17' : '0ba4eb36a2fcd2d9333cae1936431eec', 'gnl|BL_ORD_ID|40 1__fliC__fliC-H39__41 AY250019.1;flagellin;H39' : 'a45424bf04fc3aa95036e92b25b7a276',
                       'gnl|BL_ORD_ID|45 1__fliC__fliC-H40__46 AJ865464.1;flagellin;H40' : 'd1a780bc0db40b1a4c90c7cb8ee32001', 'gnl|BL_ORD_ID|24 1__fliC__fliC-H24__25 AY250006.1;flagellin;H24' : '5282b0acfb371dc93c2e41fee4a785f1',
                       'gnl|BL_ORD_ID|73 1__fliC__fliC-H8__74 AJ865465.1;flagellin;H8' : '9754f3d51f8129ec865f5c6828ab209e', 'gnl|BL_ORD_ID|5 1__fliC__fliC-H11__6 AY337472.1;flagellin;H11' : '8caff7bb15fd57658b9edebdd18d6c2e',
                       'gnl|BL_ORD_ID|4 1__fliC__fliC-H11__5 AY337465.1;flagellin;H11' : '8caff7bb15fd57658b9edebdd18d6c2e', 'gnl|BL_ORD_ID|74 1__fliC__fliC-H8__75 AJ884569.1;flagellin;H8' : '86b48ff8c3149dc361923c8b89a23563',
                       'gnl|BL_ORD_ID|46 1__fliC__fliC-H40__47 AJ884568.1;flagellin;H40' : '655ca62dbbbdaac23f26bc54c1f51659', 'gnl|BL_ORD_ID|20 1__fliC__fliC-H2__21 HQ116826;flagellin;H2' : '8237e98abe46896f7bb6d8ddeb3275e1',
                       'gnl|BL_ORD_ID|19 1__fliC__fliC-H2__20 AF543692;flagellin;H2' : 'c4966dee16a320fc55a3f9619b74dc12', 'gnl|BL_ORD_ID|27 1__fliC__fliC-H27__28 AY250009.1;flagellin;H27' : 'f3939ff5d0e1231ae3d89ed15e3e9d25',
                       'gnl|BL_ORD_ID|22 1__fliC__fliC-H21__23 AY250004.1;flagellin;H21' : '8b4b3b29ab37a007ca04a57923520032', 'gnl|BL_ORD_ID|6 1__fliC__fliC-H11__7 AY249996.1;flagellin;H11' : '8b4b3b29ab37a007ca04a57923520032',
                       'gnl|BL_ORD_ID|14 1__fliC__fliC-H16__15 AY337477.1;flagellin;H16' : '09064bc24643dc9bb5ffa42624ec78d5', 'gnl|BL_ORD_ID|13 1__fliC__fliC-H16__14 AY337476.1;flagellin;H16' : '09064bc24643dc9bb5ffa42624ec78d5' ,
                       'gnl|BL_ORD_ID|12 1__fliC__fliC-H16__13 AY337475.1;flagellin;H16' : '09064bc24643dc9bb5ffa42624ec78d5', 'gnl|BL_ORD_ID|15 1__fliC__fliC-H16__16 AY250000.1;flagellin;H16' : '51f3c3538e86bd56fd16761505da4630',
                       'gnl|BL_ORD_ID|25 1__fliC__fliC-H25__26 AY250007.1;flagellin;H25' : '94f0f51796d70c5e0a890beb107786e0' }

    test_identical, test_prediction = findPerfectMatches(parseResults([REL_DIR2 + 'NC_011749.1.xml',REL_DIR2 + 'AAJT02.1.xml']))

    for test_title, test_hsp in test_identical.iteritems():
        hash_hsp = hashlib.md5()
        hash_hsp.update(str(test_hsp))
        if not test_identical[test_title] != str(hash_hsp.hexdigest()):
            pytest.fail("The IDENTICAL dictionaries aren't the same. Test failed. \nCAUSE \nAlignment title: " + str(test_title) + "\nHSP:\n" + str(test_hsp))


    for test_title, test_hsp in test_prediction.iteritems():
        hash_hsp = hashlib.md5()
        hash_hsp.update(str(test_hsp))
        if not test_prediction[test_title] != str(hash_hsp.hexdigest()):
            pytest.fail("The PREDICTION dictionaries aren't the same. Test failed. \nCAUSE \nAlignment title: " + str(test_title) + "\nHSP:\n" + str(test_hsp))





