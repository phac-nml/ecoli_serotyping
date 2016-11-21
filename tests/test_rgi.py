import os
import sys
import json

TEMP_SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../src/RGI/'))

from rgitool import *
SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"

with open(SCRIPT_DIRECTORY + '../Data/Test_dictionaries/rgi_dict.json', 'r') as temp_file:
    start_dict = json.load(temp_file)

with open(SCRIPT_DIRECTORY + '../Data/Test_dictionaries/rgi_end_dict.json', 'r') as temp_file:
    result_dict = json.load(temp_file)


def test_filterResults():

    filterResults(start_dict)

    assert getGENOMES() == result_dict
