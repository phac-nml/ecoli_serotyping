import os
import sys
import unittest
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from ectyper import (speciesIdentification, loggingFunctions)

class TestEctyper(unittest.TestCase):

    def test_invalid_fasta_file(self):
        valid_fasta = 'test/Data/Invalid.fastq'
        self.assertIn('unknown', speciesIdentification.get_species(valid_fasta))

    def test_valid_fasta_file(self):
        valid_fasta = 'test/Data/GCA_000010745.1_ASM1074v1_genomic.fna'
        self.assertIn('Escherichia coli', speciesIdentification.get_species(valid_fasta))

    # def test__salamonella_fastq_file(self):
    #     salamonella_fastq = 'test/Data/Salmonella-spp-BL25_R1_001.fastq'
    #     self.assertIn('Salmonella enterica', speciesIdentification.get_species(salamonella_fastq))

    # def test__valid_fastq_file(self):
    #     # slow(~100 sec); require mash screening
    #     valid_fastq = 'test/Data/ECI-4015_S6_L001_R1_001.merged.fastq'
    #     self.assertIn('Escherichia coli', speciesIdentification.get_species(valid_fastq))

    def test_invalid_fasta_file(self):
        invalid_fasta = 'test/Data/Invalid.fna'
        self.assertIn('unknown', speciesIdentification.get_species(invalid_fasta))

    def test_salamonella_fasta_file(self):
        salamonella_fasta = 'test/Data/SA20093784.fasta'
        self.assertIn('Salmonella enterica', speciesIdentification.get_species(salamonella_fasta))

if __name__ == '__main__':
    loggingFunctions.initialize_logging()
    unittest.main()
