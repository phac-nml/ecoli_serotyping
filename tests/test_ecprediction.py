import hashlib
import pytest
import random
import json

from ecoli_serotyping.src.Serotyper.ecprediction import *
from ecoli_serotyping.src.ecvalidatingfiles import *

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
REL_DIR = SCRIPT_DIRECTORY + '../temp/'

with open(SCRIPT_DIRECTORY + '../Data/Test_dictionaries/ecprediction_dict.json') as f:
    expected_dict = json.load(f)


initializeDB()

runBlastQuery(sorted([
    SCRIPT_DIRECTORY + '../Data/Testing_Data/Reference_Genomes/NC_011749.1.fasta'
    , SCRIPT_DIRECTORY + '../Data/Testing_Data/Reference_Genomes/AAJT02.1.fsa_nt']), 'ECTyperDB')

def test_filterPredictions():

     test_prediction = parseResults(REL_DIR + 'xml/combined_genomes.xml')
     i = 0

     while i<4:
        id_num = random.random() * 100
        length_num = random.random() * 100
        exp_dict = {}
        length_exp_dict = 0

        for expected_genome, expected_value in expected_dict.iteritems():
            if isinstance(expected_value, dict):
                expected_temp = {}
                for expected_title, expected_list in expected_value.iteritems():
                    expected_id_num = expected_list[2] * 100
                    expected_length_num = expected_list[1] * 100
                    if expected_length_num >= length_num:
                        if expected_id_num >=id_num:
                            expected_temp[expected_title] = [expected_list[0], expected_list[3]]
                            exp_dict[expected_genome]= expected_temp
                            length_exp_dict +=1
                        else:
                            expected_temp[expected_title] = [expected_list[0], 'NA']
                            exp_dict[expected_genome] = expected_temp
                            length_exp_dict +=1
                    else:
                        expected_temp[expected_title] = [expected_list[0], 'NA']
                        exp_dict[expected_genome] = expected_temp
                        length_exp_dict +=1
            else:
                exp_dict[expected_genome] = 'NA'
                length_exp_dict +=1

        test_pred = filterPredictions(test_prediction, id_num, length_num)

        length_test_pred = 0

        for test_genome, test_value in test_pred.iteritems():
            if test_value != 'NA':
                for test_title, test_hsp in test_value.iteritems():
                    hash_hsp = hashlib.md5()
                    hash_hsp.update(str(test_hsp[0]))
                    if exp_dict[test_genome][test_title][0] != str(hash_hsp.hexdigest()):
                        pytest.fail("The dictionaries aren't the same. Test failed. \nCAUSE \nAlignment title: " + str(test_title)  + '\nExpected: ' + str(hash_hsp.hexdigest()) +  '\nFound: ' +
                                         str(exp_dict[test_genome][test_title][0]))
                    if exp_dict[test_genome][test_title][1] != test_hsp[1]:
                        pytest.fail("The dictionaries aren't the same. Test failed. \nCAUSE \nAlignment title: " + str(test_title) + "\nExpected length:\n" + str(test_hsp[1]) + '\nFound length:' +
                                    str(exp_dict[test_genome][test_title][1]))

                    length_test_pred += 1
            else:
                if exp_dict[test_genome] != 'NA':
                    pytest.fail("The dictionaries aren't the same. Test failed. \nCAUSE \nGenome: \n" + test_genome)
                    length_test_pred +=1

        if length_exp_dict != length_test_pred:
            pytest.fail("The resulting dictionary is not the same length as the expected dictionary.")
        i+=1


def test_getProductPercentage():

    test_prediction = parseResults(REL_DIR + 'xml/combined_genomes.xml')
    exp_dict = {}

    for expected_genome, expected_value in expected_dict.iteritems():
        expected_temp = {}
        for expected_title, expected_list in expected_value.iteritems():
            expected_id_num = expected_list[2]
            expected_length_num = expected_list[1]
            exp_product = expected_id_num * expected_length_num
            expected_temp[expected_title] = "{0:.4f}".format(exp_product)
        exp_dict[expected_genome]  = expected_temp

    for test_genome, test_value in test_prediction.iteritems():
        if test_value != 'NA':
            for test_title, test_info in test_value.iteritems():
                test_length = test_info.keys()[0]
                test_hsp = test_info[test_length]
                test_product = getProductPercentage(test_length, test_hsp)
                test_product = "{0:.4f}".format(test_product)

                assert getProductPercentage(1, test_hsp) == 0

                assert getProductPercentage(0, test_hsp) == 0

                assert test_product == exp_dict[test_genome][test_title]

def test_sortMatches():

    test_sortmatches = sortMatches(filterPredictions(parseResults(REL_DIR + 'xml/combined_genomes.xml'), 0, 0))

    for test_genome, test_typedict in test_sortmatches.iteritems():
        for test_type, test_infolist in test_typedict.iteritems():
            for test_info in test_infolist:
                test_title = test_info.get('title')
                assert expected_dict[test_genome][test_title][4] == test_type

    assert sortMatches({'NC_011749': 'NA'}) == {'NC_011749': 'NA'}


def test_searchType():

    for test_genome, test_value in expected_dict.iteritems():
        for test_title, test_info in test_value.iteritems():
            assert searchType(test_title, str(test_info[5])[0]) == test_info[5]

    assert searchType('fli A456', 'O') == 'none'
    assert searchType(' ', 'O') == 'none'
    assert searchType('fliC-H15', ' ') == 'none'


def test_findTopMatches():

    exp_top_match1 = {'AAJT0200':{
        'gnl|BL_ORD_ID|377 9__wzy__wzy-O148__378 DQ167407.1;O antigen polyermase;O148' : ['c10c9586ec7845c902fbbd55f7df991f', 1.0, 1.0, 1143, 'wzy', 'O148'],
        #'gnl|BL_ORD_ID|28 1__fliC__fliC-H28__29 AY250010.1;flagellin;H28' : ['6dfdab8081e7d834f1b3243985472b2d', 0.987931034483, 0.987931034483, 1740, 'fliC', 'H28'],
        'gnl|BL_ORD_ID|622 fliC_105_AAJT02000052_H28': ['4ef728baef07d5a2855781ad76654218', 1.0, 1.0, 1740, 'fliC', 'H28']
        }}

    exp_top_match_2 = {'GENOME_1': 'NA', 'GENOME_2':  {'wzy': [{'title': 'gnl|BL_ORD_ID|377 9__wzy__wzy-O148__378 DQ167407.1;O antigen polyermase;O148', 'hsp': 'example_hsp', 'length': '1166', 'perc': 0.964}]},
                       'GENOME_3': {'wzx':[{'title': 'gnl|BL_ORD_ID|377 9__wzx__wzx-O148__378 DQ167407.1;O antigen polyermase;O148', 'hsp': 'example_hsp', 'length': 'NA', 'perc': 0}]}}

    test_prediction = (filterPredictions(parseResults(REL_DIR + 'xml/combined_genomes.xml'), 0, 0))
    test_top_matches = findTopMatches(sortMatches(test_prediction))

    for test_genome, test_value in test_top_matches.iteritems():
        test_otype = test_value.get('otype')
        test_htype = test_value.get('htype')
        otype_title = test_otype['title']
        htype_title = test_htype['title']

        assert exp_top_match1[test_genome][otype_title][5] == 'O148'
        assert exp_top_match1[test_genome][htype_title][5] == 'H28'

    assert findTopMatches(exp_top_match_2) == {
        'GENOME_1': {'otype': 'NA', 'htype': 'NM', 'prediction_strength': 'NA'},
        'GENOME_2': {'otype': {'title': 'gnl|BL_ORD_ID|377 9__wzy__wzy-O148__378 DQ167407.1;O antigen polyermase;O148', 'hsp': 'example_hsp', 'length': '1166', 'perc': 0.964}, 'htype': 'NM', 'prediction_strength': 'Top match'},
        'GENOME_3': {'otype': 'NA', 'htype': 'NM', 'prediction_strength': 'Top match'}
    }