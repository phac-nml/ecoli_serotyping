import os
import sys
import json

TEMP_SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../src/Tools_Controller/'))
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../src/'))

from tools_controller import *

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"

with open(SCRIPT_DIRECTORY + '../Data/Test_dictionaries/rgi_filter_results_end_dict.json') as f:
    amr_dict = json.load(f)

with open(SCRIPT_DIRECTORY + '../Data/Test_dictionaries/tools_vf_results_dict.json') as f:
    vf_dict = json.load(f)

with open(SCRIPT_DIRECTORY + '../Data/Test_dictionaries/tools_serotyper_results_dict.json') as f:
    serotyper_dict = json.load(f)

with open(SCRIPT_DIRECTORY + '../Data/Test_dictionaries/tools_serotyper_results_full_info_dict.json') as f:
    serotyper_dict_full = json.load(f)


def test_mergeResults():

    with open(SCRIPT_DIRECTORY+ '../Data/Test_dictionaries/tools_serotyper+vf_with_verbose_dict.json') as f:
        temp_dict = json.load(f)
    assert mergeResults(serotyper_dict_full, vf_dict, {}) == temp_dict

    with open(SCRIPT_DIRECTORY+ '../Data/Test_dictionaries/tools_vf+amr_dict.json') as f:
        temp_dict = json.load(f)
    assert mergeResults({}, vf_dict, amr_dict) == temp_dict

    with open(SCRIPT_DIRECTORY+ '../Data/Test_dictionaries/tools_vf_dict.json') as f:
        temp_dict = json.load(f)
    assert mergeResults({}, vf_dict, {}) == temp_dict

    with open(SCRIPT_DIRECTORY+ '../Data/Test_dictionaries/tools_all_dict.json') as f:
        temp_dict = json.load(f)
    assert mergeResults(serotyper_dict, vf_dict, amr_dict) == temp_dict