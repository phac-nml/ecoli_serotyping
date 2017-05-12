
import src.serotypePrediction
from argparse import Namespace
import json
from definitions import ROOT_DIR


CONST_CSV = False
CONST_INPUT = 'filename'
CONST_MINGENOMES = 1
CONST_PERIDENT = 90
CONST_PERLEN = 90
CONST_SER = True
CONST_VIR = True

with open(ROOT_DIR + '/Data/test_dictionaries/serotype_test_data.json') as f:
    json_output = json.load(f)


args = Namespace(csv = CONST_CSV, input = CONST_INPUT, minimumGenomes = CONST_MINGENOMES, percentIdentity =CONST_PERIDENT, percentLength = CONST_PERLEN, serotyper = CONST_SER,  virulenceFactors = CONST_VIR)

def test_parse_serotype() :

    for entry in json_output:
        for record in json_output[entry]['serotype']:
            comparison = json_output[entry]['serotype'][record]
            result = src.serotypePrediction.parse_serotype(json_output[entry]['serotype'][record]['blast_record'], args)['serotype']
            assert(record in result)
            result = result[record]

            assert(comparison['antigen'] == result['antigen'])
            result = result['blast_record']
            comparison = comparison['blast_record']

            for item in result:
                assert(result[item] == comparison[item])


def test_predict_serotype() :

    for entry in json_output:
        comparison = {}
        comparison[entry] = json_output[entry]
        result = src.serotypePrediction.predict_serotype(comparison)


        for item in comparison:
            assert(result[entry]['otype']['strength'] == comparison[entry]['otype']['strength'])
            assert (result[entry]['otype']['ant_number'] == comparison[entry]['otype']['ant_number'])

            assert(result[entry]['htype']['strength'] == comparison[entry]['htype']['strength'])
            assert (result[entry]['htype']['ant_number'] == comparison[entry]['htype']['ant_number'])

