import unittest
import argparse
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

import src.ectyper as ectyper

class TestEctyper(unittest.TestCase):

    def test_reads_option(self):
        incorrect_file_format = ['-i', 'Data/test_reads/madeup.fq','-r']
        correct_file_format = ['-i', 'Data/test_reads/ECI-4015_S6_L001_R1_001.merged.fastq','-r']

        with self.assertRaises(SystemExit):
            sys.argv[1:] = incorrect_file_format
            ectyper.run_program()

        # sys.argv[1:] = correct_file_format
        # ectyper.run_program()

if __name__ == '__main__':
    unittest.main()
    