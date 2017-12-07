import unittest, sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from ectyper.genomeFunctions import get_files_as_list, get_valid_format

class TestGenomeFunction(unittest.TestCase):

    def test_get_files_as_list(self):
        input = 'test/Data/test_dir/sample.fasta,test/Data/test_dir/sample.fasta.tar,test/Data/test_dir/sample2.fasta,test/Data/test_dir/test_junk.txt'
        output = ''
        self.assertEqual(len(get_files_as_list(input)), 3)

    def test_get_valid_format(self):
        self.assertEqual(None, get_valid_format('123'))
        self.assertEqual(None, get_valid_format('test/Data/test_dir/test_junk.txt'))
        # self.assertEqual(None, get_valid_format('test/Data/test_dir/sample.fasta.tar'))
if __name__ == '__main__':
    unittest.main()
