import os
import sys
import unittest
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from ectyper import (speciesIdentification, loggingFunctions, definitions)

class TestSpeciesId(unittest.TestCase):
    
    def test_salamonella_fastq_file(self):
        if not os.path.isfile(definitions.REFSEQ_SKETCH):
            print("No seqref file. Skip this test.")
            return
        print("Testing non-ecoli fastq")
        salamonella_fastq = 'test/Data/Salmonella.fastq'
        self.assertIn('Salmonella enterica', speciesIdentification.get_species(salamonella_fastq))

    def test_ecoli_fastq_file(self):
        if not os.path.isfile(definitions.REFSEQ_SKETCH):
            print("No seqref file. Skip this test.")
            return
        print("Testing ecoli fastq")
        valid_fastq = 'test/Data/Escherichia.fastq'
        self.assertIn('Escherichia coli', speciesIdentification.get_species(valid_fastq))

    def test_ecoli_fasta_file(self):
        if not os.path.isfile(definitions.REFSEQ_SKETCH):
            print("No seqref file. Skip this test.")
            return
        print("Testing ecoli fasta")
        if not os.path.isfile(definitions.REFSEQ_SKETCH):
            return
        valid_fasta = 'test/Data/Escherichia.fna'
        self.assertIn('Escherichia coli', speciesIdentification.get_species(valid_fasta))

    def test_different_species_fasta_file(self):
        if not os.path.isfile(definitions.REFSEQ_SKETCH):
            print("No seqref file. Skip this test.")
            return
        print("Testing non-ecoli fasta")
        salamonella_fasta = 'test/Data/Salmonella.fasta'
        streptococcus_fasta = 'test/Data/Streptococcus.fasta'
        straphylococcus_fasta = 'test/Data/Straphylococcus.fasta'
        yersinia_fasta = 'test/Data/Yersinia.fasta'
        listeria_fasta = 'test/Data/Listeria.fasta'
        campylobacter_fasta = 'test/Data/Campylobacter.fasta'
        self.assertIn('Salmonella', speciesIdentification.get_species(salamonella_fasta))
        self.assertIn('Streptococcus', speciesIdentification.get_species(streptococcus_fasta))
        self.assertIn('Staphylococcus', speciesIdentification.get_species(straphylococcus_fasta))
        self.assertIn('Yersinia', speciesIdentification.get_species(yersinia_fasta))
        self.assertIn('Listeria', speciesIdentification.get_species(listeria_fasta))
        self.assertIn('Campylobacter', speciesIdentification.get_species(campylobacter_fasta))

if __name__ == '__main__':
    loggingFunctions.initialize_logging()
    unittest.main()
