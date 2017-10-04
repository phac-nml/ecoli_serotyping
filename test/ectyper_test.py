import unittest
import argparse
import sys
import os

from ectyper import ectyper

class TestEctyper(unittest.TestCase):

    def test_fastq_file(self):
        salamonella_fastq = ['-i', 'test/Data/Salmonella.fastq','-r']
        valid_fastq = ['-i', 'test/Data/Ecoli_O22H8.fastq','-r']
        invalid_fastq = ['-i', 'test/Data/invalid.fastq','-r']

        with self.assertRaises(SystemExit):
            sys.argv[1:] = invalid_fastq
            ectyper.run_program()

        with self.assertRaises(SystemExit):
            sys.argv[1:] = salamonella_fastq
            ectyper.run_program()

        sys.argv[1:] = valid_fastq
        ectyper.run_program()

if __name__ == '__main__':
    unittest.main()
    