from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import open
from builtins import str
from future import standard_library
standard_library.install_aliases()
import argparse
import hashlib
import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from ectyper import ectyper




def set_input(input, percent_iden=None, output=None):
    args = ['-i', input]
    if percent_iden:
        args += ['-d', str(percent_iden)]
    if output:
        args += ['-o', output]
    sys.argv[1:] = args


def get_md5(file):
    with open(file, 'rb') as handler:
        return hashlib.md5(handler.read()).hexdigest()


class TestEctyper(unittest.TestCase):
    # Print all differences if there is difference
    maxDiff = None

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
        input = 'test/Data/Escherichia.fastq'
        output_dir = 'test_valid_fastq_file'
        expected_checksum = '95d5cc69090fc686d3f5c5ba6602f374'
        set_input(input, output=output_dir)
        ectyper.run_program()
        self.assertEqual(
            get_md5(os.path.join('output', output_dir, 'output.csv')),
            expected_checksum
        )

    def test_folder_input(self):
        input = 'test/Data/test_dir'
        output_dir = 'test_folder_input'
        set_input(input, output=output_dir)
        expected_checksum = 'bd9354495c264bcfb79c5581b9e7c953'
        ectyper.run_program()
        self.assertEqual(
            get_md5(os.path.join('output', output_dir, 'output.csv')),
            expected_checksum
        )

    def test_list_input(self):
        input = 'test/Data/test_dir/sample.fasta,test/Data/test_dir/sample.fasta'
        output_dir = 'test_list_input'
        set_input(input, output=output_dir)
        expected_checksum = '3ad690bf9f3ea87dd0200600f72fa618'
        ectyper.run_program()
        self.assertEqual(
            get_md5(os.path.join('output', output_dir, 'output.csv')),
            expected_checksum
        )

    def test_create_tmp_files(self):
        # Test relative path
        expected_dict = {'assemble_temp_dir': 'test/temp/assemblies',
                         'fasta_temp_dir': 'test/temp/fastas',
                         'output_dir': os.path.abspath('output') + '/',
                         'output_file': os.path.abspath('output/output.csv')}
        self.assertDictEqual(expected_dict, ectyper.create_tmp_files(
            'test/temp', output_dir=''))
        expected_dict = {'assemble_temp_dir': 'test/temp/assemblies',
                         'fasta_temp_dir': 'test/temp/fastas',
                         'output_dir': os.path.abspath('output/test'),
                         'output_file': os.path.abspath('output/test/output.csv')}
        self.assertDictEqual(expected_dict, ectyper.create_tmp_files(
            'test/temp', output_dir='test'))


if __name__ == '__main__':
    unittest.main()
