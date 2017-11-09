import argparse
import os
import sys
import unittest
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from ectyper import ectyper

def set_input(input, allow_species_id=True, percent_iden=97):
    args = ['-i', input, '-d', str(percent_iden)]
    if allow_species_id:
        args.append('-s')
    sys.argv[1:] = args

class TestEctyper(unittest.TestCase):
    '''
    Integration tests (slow)
    '''

    def test__valid_fastq_file(self):
        valid_fastq = 'test/Data/Escherichia.fastq'
        set_input(valid_fastq)
        ectyper.run_program()

    # def test__valid_fasta_file(self):
    #     valid_fasta = 'test/Data/Escherichia.fna'
    #     set_input(valid_fasta)
    #     ectyper.run_program()
if __name__ == '__main__':
    unittest.main()
