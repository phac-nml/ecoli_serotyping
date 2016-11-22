import os
import sys
import json
import pytest

TEMP_SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../src/Virulence_Factors/'))
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../src/'))

from virulencefactors import *
from sharedmethods import *

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"

with open(SCRIPT_DIRECTORY + '../Data/Test_dictionaries/virulencefactors_dict_expected.json') as f:
    expectedDict = json.load(f)

with open(SCRIPT_DIRECTORY + '../Data/Test_dictionaries/virulencefactors_dict_min2.json') as f:
    resultDict1 = json.load(f)

with open(SCRIPT_DIRECTORY + '../Data/Test_dictionaries/virulencefactors_dict_min6.json') as f:
    resultDict2 = json.load(f)


def test_filterVFs():

    assert filterGenes(expectedDict, 1) == expectedDict

    assert filterGenes(expectedDict, 2) == resultDict1

    assert filterGenes(expectedDict, 6) == resultDict2

    assert filterGenes(expectedDict, 8) == {'AAJT0200': {},
                                          'NC_002127': {},
                                          'NC_002128': {},
                                          'NC_002695': {},
                                          'NC_011739': {},
                                          'NC_011749': {},
                                          'NC_011750': {},
                                          'NC_011751': {}}
