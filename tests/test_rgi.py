import os
import sys
import json

TEMP_SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../src/RGI/'))
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../src/'))

from sharedmethods import createDirs
from rgitool import *

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"

createDirs()

with open(SCRIPT_DIRECTORY + '../Data/Test_dictionaries/rgi_filter_results_start_dict.json', 'r') as temp_file:
    start_dict = json.load(temp_file)

with open(SCRIPT_DIRECTORY + '../Data/Test_dictionaries/rgi_getresults_dict.json', 'r') as temp_file:
    result_dict1 = json.load(temp_file)

with open(SCRIPT_DIRECTORY + '../Data/Test_dictionaries/rgi_filter_results_end_dict.json', 'r') as temp_file:
    result_dict2 = json.load(temp_file)



RGIpath = ''
for relpath in sys.path:
    if 'release-rgi-v3.1.1-58cad6a3b443abb290cf3df438fe558bc5bfec39' in str(relpath):
        RGIpath = relpath + '/'
        break

genomesList = [SCRIPT_DIRECTORY + '../Data/Testing_Data/Reference_Genomes/AAJT02.1.fsa_nt', SCRIPT_DIRECTORY + '../Data/Testing_Data/Reference_Genomes/NC_011749.1.fasta' ]


def test_getResults():
    clearGlobalDicts()
    checkFiles(genomesList)

    genomes = getGENOMES()
    genomeNames = getGENOMENAMES()
    clearGlobalDicts()

    setGlobalDicts(genomes, genomeNames)

    getResults(genomesList, RGIpath)

    assert getGENOMESDict() == result_dict1


def test_filterResults():

    getGeneDict(start_dict)

    assert getGENOMESDict() == result_dict2
