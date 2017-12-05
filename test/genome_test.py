import unittest, sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from ectyper.genomeFunctions import get_files_as_list, get_valid_format

class TestGenomeFunction(unittest.TestCase):

    def test_get_files_as_list(self):
        input = '/home/sam/Projects/galaxy/database/files/000/dataset_85.dat,/home/sam/Projects/galaxy/database/files/000/dataset_84.dat,/home/sam/Projects/galaxy/database/files/000/dataset_83.dat,/home/sam/Projects/galaxy/database/files/000/dataset_82.dat,/home/sam/Projects/galaxy/database/files/000/dataset_81.dat,/home/sam/Projects/galaxy/database/files/000/dataset_80.dat,/home/sam/Projects/galaxy/database/files/000/dataset_79.dat,/home/sam/Projects/galaxy/database/files/000/dataset_78.dat'
        output = ''
        self.assertEqual(len(get_files_as_list(input)), 8)

    def test_get_valid_format(self):
        get_valid_format('123')
if __name__ == '__main__':
    unittest.main()
