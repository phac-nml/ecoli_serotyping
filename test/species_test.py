from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from ectyper import definitions, speciesIdentification


class TestSpeciesId(unittest.TestCase):

    def test_ecoli_fasta_file(self):
        escherichia_fasta = 'test/Data/Escherichia.fna'
        salamonella_fasta = 'test/Data/Salmonella.fasta'
        self.assertTrue(
            speciesIdentification.is_ecoli_genome(escherichia_fasta))
        self.assertFalse(
            speciesIdentification.is_ecoli_genome(salamonella_fasta))

    def test_different_species_fasta_file(self):
        if not os.path.isfile(definitions.REFSEQ_SKETCH):
            print("No seqref file. Skip this test.")
            return False
        escherichia_fasta = 'test/Data/Escherichia.fna'
        salamonella_fasta = 'test/Data/Salmonella.fasta'
        streptococcus_fasta = 'test/Data/Streptococcus.fasta'
        straphylococcus_fasta = 'test/Data/Straphylococcus.fasta'
        yersinia_fasta = 'test/Data/Yersinia.fasta'
        listeria_fasta = 'test/Data/Listeria.fasta'
        campylobacter_fasta = 'test/Data/Campylobacter.fasta'
        self.assertIn('Escherichia coli',
                      speciesIdentification.get_species(escherichia_fasta))
        self.assertIn(
            'Salmonella', speciesIdentification.get_species(salamonella_fasta))
        self.assertIn('Streptococcus',
                      speciesIdentification.get_species(streptococcus_fasta))
        self.assertIn('Staphylococcus', speciesIdentification.get_species(
            straphylococcus_fasta))
        self.assertIn(
            'Yersinia', speciesIdentification.get_species(yersinia_fasta))
        self.assertIn(
            'Listeria', speciesIdentification.get_species(listeria_fasta))
        self.assertIn('Campylobacter',
                      speciesIdentification.get_species(campylobacter_fasta))


if __name__ == '__main__':
    unittest.main()
