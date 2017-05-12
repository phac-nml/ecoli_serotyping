import src.genomeFunctions
from definitions import ROOT_DIR



def test_files_as_list() :

    test_list = [ROOT_DIR + '/Data/test_genomes/GCA_000010745.1_ASM1074v1_genomic.fna', ROOT_DIR + '/Data/test_genomes/GCA_000017745.1_ASM1774v1_genomic.fna', ROOT_DIR + '/Data/test_genomes/GCA_000025165.1_ASM2516v1_genomic.fna', ROOT_DIR + '/Data/test_genomes/GCA_000026545.1_ASM2654v1_genomic.fna', ROOT_DIR + '/Data/test_genomes/GCA_000091005.1_ASM9100v1_genomic.fna']

    assert(src.genomeFunctions.get_files_as_list('') == [])
    assert(src.genomeFunctions.get_files_as_list(ROOT_DIR + '/Data/test_genomes/empty') == [])
    assert (src.genomeFunctions.get_files_as_list(ROOT_DIR + '/Data/test_genomes') == test_list)
    assert (src.genomeFunctions.get_files_as_list(ROOT_DIR + '/Data/test_genomes/GCA_000010745.1_ASM1074v1_genomic.fna') ==test_list[0:1])

