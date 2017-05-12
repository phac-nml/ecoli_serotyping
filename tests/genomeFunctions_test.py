import src.genomeFunctions
from definitions import ROOT_DIR

TEST_LIST = [ROOT_DIR + '/Data/test_genomes/GCA_000010745.1_ASM1074v1_genomic.fna',
             ROOT_DIR + '/Data/test_genomes/GCA_000017745.1_ASM1774v1_genomic.fna',
             ROOT_DIR + '/Data/test_genomes/GCA_000025165.1_ASM2516v1_genomic.fna',
             ROOT_DIR + '/Data/test_genomes/GCA_000026545.1_ASM2654v1_genomic.fna',
             ROOT_DIR + '/Data/test_genomes/GCA_000091005.1_ASM9100v1_genomic.fna']


TEST_LIST2 = [ROOT_DIR + '/Data/test_genomes/bad_file.fna']

TEST_LIST_NAMES = ['AP010958.1', 'CP000800.1', 'CP001846.1', 'FM180568.1', 'AP010953.1']

def test_files_as_list() :
    assert(src.genomeFunctions.get_files_as_list('') == [])
    assert(src.genomeFunctions.get_files_as_list(ROOT_DIR + '/Data/test_genomes/empty') == [])
    assert (src.genomeFunctions.get_files_as_list(ROOT_DIR + '/Data/test_genomes') == TEST_LIST)
    assert (src.genomeFunctions.get_files_as_list(ROOT_DIR + '/Data/test_genomes/GCA_000010745.1_ASM1074v1_genomic.fna') == TEST_LIST[0:1])



def test_validate_fasta_files():

    assert(src.genomeFunctions.validate_fasta_files(TEST_LIST) == TEST_LIST)
    assert(src.genomeFunctions.validate_fasta_files([])) == []
    assert src.genomeFunctions.validate_fasta_files(TEST_LIST2) == []


#def test_get_genome_names_from_files():

   #assert src.genomeFunctions.get_genome_names_from_files(TEST_LIST) == TEST_LIST_NAMES, TEST_LIST


#est_get_genome_names_from_files()