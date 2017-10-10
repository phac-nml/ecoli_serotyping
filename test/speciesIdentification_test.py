import os
import sys
import unittest
import logging
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from ectyper import (speciesIdentification, loggingFunctions, definitions)

LOG = logging.getLogger(__name__)

class TestEctyper(unittest.TestCase):

    def test_invalid_fastq_file(self):
        LOG.info("Testing invalid fastq")
        valid_fasta = 'test/Data/Invalid.fastq'
        self.assertIn('unknown', speciesIdentification.get_species(valid_fasta))
    
    def test__salamonella_fastq_file(self):
        LOG.info("Testing non-ecoli fastq")
        salamonella_fastq = 'test/Data/Salmonella-spp-BL25_R1_001.fastq'
        self.assertIn('Salmonella enterica', speciesIdentification.get_species(salamonella_fastq))

    def test__valid_fastq_file(self):
        # slow(~100 sec); require mash screening
        LOG.info("Testing ecoli fastq")
        valid_fastq = 'test/Data/ECI-4015_S6_L001_R1_001.merged.fastq'
        self.assertIn('Escherichia coli', speciesIdentification.get_species(valid_fastq))

    def test_valid_fasta_file(self):
        LOG.info("Testing ecoli fasta")
        if not os.path.isfile(definitions.REFSEQ_SKETCH):
            return
        valid_fasta = 'test/Data/GCA_000010745.1_ASM1074v1_genomic.fna'
        self.assertIn('Escherichia coli', speciesIdentification.get_species(valid_fasta))


    def test_invalid_fasta_file(self):
        LOG.info("Testing invalid fasta")
        invalid_fasta = 'test/Data/Invalid.fna'
        self.assertIn('unknown', speciesIdentification.get_species(invalid_fasta))

    def test_different_species_fasta_file(self):
        LOG.info("Testing non-ecoli fasta")
        salamonella_fasta = 'test/Data/SA20093784.fasta'
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
    if not os.path.isfile(definitions.REFSEQ_SKETCH):
        LOG.info("No seqref file. Skip this test suite.")
        exit
    unittest.main()
