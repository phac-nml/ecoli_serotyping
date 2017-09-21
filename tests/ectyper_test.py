import unittest
import argparse
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

import src.ectyper as ectyper

class TestEctyper(unittest.TestCase):

    def test_reads_option(self):
        incorrect_file_format = ['-i', 'Data/test_reads/lambda_virus.fa','-r']
        correct_file_format = ['-i', 'Data/test_reads/reads_1.fq','-r','-s']

        with self.assertRaises(RuntimeError):
            sys.argv[1:] = incorrect_file_format
            ectyper.run_program()

        sys.argv[1:] = correct_file_format
        ectyper.run_program()

if __name__ == '__main__':
    unittest.main()
    