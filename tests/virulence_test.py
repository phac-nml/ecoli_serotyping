import src.virulencePrediction
import pandas as pd
import os
import sys
import json
from argparse import Namespace
from definitions import ROOT_DIR


CONST_CSV = False
CONST_INPUT = 'filename'
CONST_MINGENOMES = 1
CONST_PERIDENT = 80
CONST_PERLEN = 90
CONST_SER = True
CONST_VIR = True

#need to make tests for edge case (no virulence factors)



output = {}

#used to generate data
"""
for entry in test_data:
    genome = {}
    for factor in test_data[entry]:
        output[entry] = ''
        for spec in test_data[entry][factor]:
            genome[spec] = 1
            output[entry] = genome
"""

args = Namespace(csv = CONST_CSV, input = CONST_INPUT, minimumGenomes = CONST_MINGENOMES, percentIdentity =CONST_PERIDENT, percentLength = CONST_PERLEN, serotyper = CONST_SER,  virulenceFactors = CONST_VIR)

with open(ROOT_DIR + '/Data/test_dictionaries/virulence_factor_test_data.json') as f:
    test_data = json.load(f)

with open(ROOT_DIR + '/Data/test_dictionaries/expected_virulence_test_data.json') as f:
    data = json.load(f)

def test_parse__virulence_factors() :
    for entry in data:
        for vf in test_data[entry]['vf']:
            result = src.virulencePrediction.parse_virulence_factors(test_data[entry]['vf'][vf]['blast_record'], args)
            prediction = vf
            for res in result['vf']:
                assert(res == prediction)


