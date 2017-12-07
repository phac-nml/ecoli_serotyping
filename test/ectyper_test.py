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
    def test_invalid_input(self):
        valid_fasta = 'test/Data/test_dir/sample.fasta'
        set_input(valid_fasta, percent_iden=999)
        with self.assertRaises(SystemExit) as e:
            ectyper.run_program()
        self.assertEqual(e.exception.code, 2)
    
    def test_invalid_file(self):
        invalid_file = ''
        set_input(invalid_file)
        with self.assertRaises(SystemExit) as e:
            ectyper.run_program()
        self.assertEqual(e.exception.code, 0)

    def test_valid_fastq_file(self):
        valid_fastq = 'test/Data/Escherichia.fastq'
        set_input(valid_fastq)
        ectyper.run_program()

    def test_folder_input(self):
        valid_fasta = 'test/Data/test_dir'
        set_input(valid_fasta)
        ectyper.run_program()
    
    def test_list_input(self):
        valid_fasta = 'test/Data/test_dir/sample.fasta,test/Data/test_dir/sample.fasta'
        set_input(valid_fasta)
        ectyper.run_program()
if __name__ == '__main__':
    unittest.main()
