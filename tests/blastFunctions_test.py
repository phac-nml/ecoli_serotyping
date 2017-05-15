import src.blastFunctions
from argparse import Namespace
from definitions import ROOT_DIR
import tempfile

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

def test_record_passes_cutoffs():
    assert src.blastFunctions.record_passes_cutoffs(TEST_RECORD, args) == True
    assert src.blastFunctions.record_passes_cutoffs(TEST_RECORD, args_fi) == False
    assert src.blastFunctions.record_passes_cutoffs(TEST_RECORD, args_fl) == False
    assert src.blastFunctions.record_passes_cutoffs(TEST_RECORD_FI, args) == False
    assert src.blastFunctions.record_passes_cutoffs(TEST_RECORD_FL, args) == False


def test_create_blast_db():
    assert src.blastFunctions.create_blast_db(TEST_LIST) == tempfile + '/' + 'ectyper_blastdb'
