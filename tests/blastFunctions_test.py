import src.blastFunctions
from argparse import Namespace
from definitions import ROOT_DIR
import tempfile
import definitions
import hashlib

tempfile = tempfile.gettempdir()


CONST_CSV = False
CONST_INPUT = 'filename'
CONST_MINGENOMES = 1
CONST_PERIDENT = 90
CONST_PERIDENT_F = 101
CONST_PERLEN = 90
CONST_PERLEN_F = 101
CONST_SER = True
CONST_VIR = True

args = Namespace(csv = CONST_CSV, input = CONST_INPUT, minimumGenomes = CONST_MINGENOMES, percentIdentity =CONST_PERIDENT, percentLength = CONST_PERLEN, serotyper = CONST_SER,  virulenceFactors = CONST_VIR)
#create args that have impossible percent Identity
args_fi = Namespace(csv = CONST_CSV, input = CONST_INPUT, minimumGenomes = CONST_MINGENOMES, percentIdentity =CONST_PERIDENT_F, percentLength = CONST_PERLEN, serotyper = CONST_SER,  virulenceFactors = CONST_VIR)
#create args that have imossible percent length
args_fl = Namespace(csv = CONST_CSV, input = CONST_INPUT, minimumGenomes = CONST_MINGENOMES, percentIdentity =CONST_PERIDENT, percentLength = CONST_PERLEN_F, serotyper = CONST_SER,  virulenceFactors = CONST_VIR)


#test record we know should meet cutoff
TEST_RECORD = {"sstart": "2856508", "length": "742", "qseqid": "wzy_199_AKLQ01000042_O157", "send": "2855767", "sframe": "1", "pident": "99.87", "qlen": "741", "sseqid": "AE005174.2"}
#set percent identity of match to 0
TEST_RECORD_FI = {"sstart": "2856508", "length": "742", "qseqid": "wzy_199_AKLQ01000042_O157", "send": "2855767", "sframe": "1", "pident": "0", "qlen": "741", "sseqid": "AE005174.2"}
#set length of match to 0
TEST_RECORD_FL = {"sstart": "2856508", "length": "0", "qseqid": "wzy_199_AKLQ01000042_O157", "send": "2855767", "sframe": "1", "pident": "99.87", "qlen": "741", "sseqid": "AE005174.2"}

TEST_LIST = [ROOT_DIR + '/Data/test_genomes/GCA_000010745.1_ASM1074v1_genomic.fna',
             ROOT_DIR + '/Data/test_genomes/GCA_000017745.1_ASM1774v1_genomic.fna',
             ROOT_DIR + '/Data/test_genomes/GCA_000025165.1_ASM2516v1_genomic.fna',
             ROOT_DIR + '/Data/test_genomes/GCA_000026545.1_ASM2654v1_genomic.fna',
             ROOT_DIR + '/Data/test_genomes/GCA_000091005.1_ASM9100v1_genomic.fna']

TEST_QUERIES = [definitions.SEROTYPE_FILE, definitions.VF_FILE]

TEST_DB = ROOT_DIR + '/Data/test_blastdb/ectyper_blastdb'
TEMP_DB = '/tmp/ectyper_blastdb'

def test_record_passes_cutoffs():
    assert src.blastFunctions.record_passes_cutoffs(TEST_RECORD, args) == True
    assert src.blastFunctions.record_passes_cutoffs(TEST_RECORD, args_fi) == False
    assert src.blastFunctions.record_passes_cutoffs(TEST_RECORD, args_fl) == False
    assert src.blastFunctions.record_passes_cutoffs(TEST_RECORD_FI, args) == False
    assert src.blastFunctions.record_passes_cutoffs(TEST_RECORD_FL, args) == False


def test_create_blast_db():
    assert src.blastFunctions.create_blast_db(TEST_LIST[0:1]) == TEMP_DB
    with open(TEMP_DB + '.nhr', 'rb') as file:
        data = file.read()
        s1 = hashlib.sha1(data).hexdigest()
        assert s1 == 'ee3f0babf9385a4305d11e98e2b199fd87703834'

    # """
    # #this hash value will always change so no test can be run for it
    # with open('/tmp/ectyper_blastdb.nin', 'rb') as file_to_check:
    #     data = file_to_check.read()
    #     s1 = hashlib.sha1(data).hexdigest()
    #     print(s1)
    # """

    with open(TEMP_DB + '.nsq', 'rb') as file:
        data = file.read()
        s1 = hashlib.sha1(data).hexdigest()
        assert s1 == 'f26f99e51f8f4f0dbb53cb47c3bdb2ed79ed8a30'


# def test_run_blast():
#     assert src.blastFunctions.run_blast(TEST_QUERIES[0], TEMP_DB) == TEMP_DB + '.output'
#     with open(ROOT_DIR + '/Data/test_blastdb/ectyper_blastdb.output', 'rb') as file:
#         data = file.read()
#         s1 = hashlib.sha1(data).hexdigest()
#         assert s1 == 'd141bec3875faa479b57a66fe643c90189ff22dd'


#
# def test_parse_blast_results():
#     list_of_dict = [{'parser': src.serotypePrediction.parse_serotype,
#                      'predictor': src.serotypePrediction.predict_serotype
#                      }, {
#                         'parser': src.virulencePrediction.parse_virulence_factors,
#                         'predictor': src.virulencePrediction.predict_virulence_factors}]
#
#
#     data_both = src.blastFunctions.parse_blast_results(args, ROOT_DIR + '/Data/test_blastdb/ectyper_blastdb.output', list_of_dict)
#
#     with open(ROOT_DIR + '/Data/test_dictionaries/blast_parse_test_data.json') as f:
#         json_output = json.load(f)
#
#     data_both = dict(data_both)
#
#     #returns an inconsistent dictionary
#     print(data_both)
#     print(json_output)
#     assert data_both == json_output
#
#     with open(ROOT_DIR + '/Data/test_dictionaries/blast_parse_ser_test_data.json') as f:
#         json_output = json.load(f)
#
#     data_ser = src.blastFunctions.parse_blast_results(args,ROOT_DIR + '/Data/test_blastdb/ectyper_blastdb.output', list_of_dict[0:1])
#     data_ser = dict(data_ser)
#     print (data_ser)
#     print(json_output)
#     assert data_ser == json_output
#
#     with open(ROOT_DIR + '/Data/test_dictionaries/blast_parse_vir_test_data.json') as f:
#         json_output = json.load(f)
#
#     data_vir = src.blastFunctions.parse_blast_results(args,ROOT_DIR + '/Data/test_blastdb/ectyper_blastdb.output',
#                                                       list_of_dict[1:2])
#     data_vir = dict(data_vir)
#
#     assert data_vir == json_output
#
#     data_none = src.blastFunctions.parse_blast_results(args,ROOT_DIR + '/Data/test_blastdb/ectyper_blastdb.output', [])
#
#     assert data_none == {}
