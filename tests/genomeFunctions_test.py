import src.genomeFunctions
from definitions import ROOT_DIR

TEST_LIST = [ROOT_DIR + '/Data/test_genomes/GCA_000010745.1_ASM1074v1_genomic.fna',
             ROOT_DIR + '/Data/test_genomes/GCA_000017745.1_ASM1774v1_genomic.fna',
             ROOT_DIR + '/Data/test_genomes/GCA_000025165.1_ASM2516v1_genomic.fna',
             ROOT_DIR + '/Data/test_genomes/GCA_000026545.1_ASM2654v1_genomic.fna',
             ROOT_DIR + '/Data/test_genomes/GCA_000091005.1_ASM9100v1_genomic.fna']


TEST_LIST2 = [ROOT_DIR + '/Data/test_genomes/bad_file.fna']

TEST_LIST_NAMES = ['AP010958.1', 'CP000800.1', 'CP001846.1', 'FM180568.1', 'AP010953.1']

TEST_HEADERS = ['AP010958.1 Escherichia coli O103:H2 str. 12009 DNA, complete genome',
                'CP000800.1 Escherichia coli E24377A, complete genome',
                'CP001846.1 Escherichia coli O55:H7 str. CB9615, complete genome',
                'FM180568.1 Escherichia coli 0127:H6 E2348/69 complete genome, strain E2348/69',
                'AP010953.1 Escherichia coli O26:H11 str. 11368 DNA, complete genome']

def test_files_as_list() :
    assert(src.genomeFunctions.get_files_as_list('') == [])
    assert (src.genomeFunctions.get_files_as_list(ROOT_DIR + '/Data/test_genomes') == TEST_LIST)
    assert (src.genomeFunctions.get_files_as_list(ROOT_DIR + '/Data/test_genomes/GCA_000010745.1_ASM1074v1_genomic.fna') == TEST_LIST[0:1])



def test_validate_fasta_files():

    assert(src.genomeFunctions.validate_fasta_files(TEST_LIST) == TEST_LIST)
    assert(src.genomeFunctions.validate_fasta_files([])) == []
    assert src.genomeFunctions.validate_fasta_files(TEST_LIST2) == []


def test_get_genome_names_from_files():
   assert src.genomeFunctions.get_genome_names_from_files(TEST_LIST) == (TEST_LIST_NAMES, TEST_LIST)
   assert src.genomeFunctions.get_genome_names_from_files(TEST_LIST[0:1]) == (TEST_LIST_NAMES[0:1], TEST_LIST[0:1])
   assert src.genomeFunctions.get_genome_names_from_files([]) == ([], [])


def test_get_fasta_header_from_file():
    assert src.genomeFunctions.get_fasta_header_from_file(ROOT_DIR + '/Data/test_genomes/bad_file.fna') == None
    assert src.genomeFunctions.get_fasta_header_from_file(TEST_LIST[0])  == TEST_HEADERS[0]

test_get_fasta_header_from_file()