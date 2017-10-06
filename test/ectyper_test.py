import argparse
import os
import sys
import unittest
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from ectyper import ectyper

def set_input(input, allow_species_id=True):
    args = ['-i', input]
    if allow_species_id:
        args.append('-s')
    sys.argv[1:] = args

class TestEctyper(unittest.TestCase):
    '''
    Integration tests (slow)
    Passed tests are commented out
    '''

    # def test__invalid_fastq_file(self):
    #     invalid_fastq = 'test/Data/Invalid.fastq'
    #     with self.assertRaises(SystemExit):
    #         set_input(invalid_fastq)
    #         ectyper.run_program()

    def test__salamonella_fastq_file(self):
        salamonella_fastq = 'test/Data/Salmonella-spp-BL25_R1_001.fastq'
        with self.assertRaises(SystemExit):
            set_input(salamonella_fastq)
            ectyper.run_program()

    # def test__valid_fastq_file(self):
    #     valid_fastq = 'test/Data/ECI-4015_S6_L001_R1_001.merged.fastq'
    #     set_input(valid_fastq)
    #     ectyper.run_program()

    # def test__invalid_fasta_file(self):
    #     invalid_fasta = 'test/Data/Invalid.fna'

    #     with self.assertRaises(SystemExit):
    #         set_input(invalid_fasta)
    #         ectyper.run_program()

    # def test__salamonella_fasta_file(self):
    #     salamonella_fasta = 'test/Data/SA20093784.fasta'
    #     with self.assertRaises(SystemExit):
    #         set_input(salamonella_fasta)
    #         ectyper.run_program()

    # def test__valid_fasta_file(self):
    #     valid_fasta = 'test/Data/GCA_000010745.1_ASM1074v1_genomic.fna'
    #     set_input(valid_fasta)
    #     ectyper.run_program()
if __name__ == '__main__':
    unittest.main()
