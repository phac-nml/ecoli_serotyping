from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from ectyper.genomeFunctions import get_files_as_list, get_valid_format


class TestGenomeFunction(unittest.TestCase):

    def test_get_files_as_list(self):
        # this function does not check for validity of file
        input = 'test/Data/test_dir/sample.fasta,test/Data/test_dir/sample.fasta.tar,test/Data/test_dir/sample2.fasta,test/Data/test_dir/test_junk.txt'
        output = ''
        self.assertEqual(len(get_files_as_list(input)), 4)

    def test_get_valid_format(self):
        # None file
        self.assertEqual(None, get_valid_format('123'))
        # Random txt
        self.assertEqual(None, get_valid_format(
            'test/Data/test_dir/test_junk.txt'))
        # renamed tar file
        self.assertEqual(None, get_valid_format(
            'test/Data/test_dir/sampletar'))
        # zip file
        self.assertEqual(None, get_valid_format(
            'test/Data/test_dir/sample.fasta.zip'))
        # actual fasta
        self.assertEqual('fasta', get_valid_format(
            'test/Data/test_dir/sample.fasta'))


if __name__ == '__main__':
    unittest.main()
