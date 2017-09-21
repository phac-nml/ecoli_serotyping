import unittest
import src.shortReads as shortReads

class TestShortReads(unittest.TestCase):

    def test_samtofasta(self):
        '''
        Test samToFasta(sam_file)
        '''
        empty_input = ''
        correct_input = 'Data/test_reads/eg1.sam'
        correct_output = 'Data/test_reads/eg1.fasta'

        with self.assertRaises(RuntimeError):
            shortReads.samToFasta(empty_input)
        self.assertTrue(shortReads.samToFasta(correct_input)==correct_output)

    def test_build_bowtie_index(self):
        '''
        Test build_bowtie_index(reference_in)
        '''
        correct_input = 'Data/test_reads/EcOH.fasta'
        correct_output = 'Data/test_reads/bowtie_index/'
        empty_input = ''
        self.assertEqual(shortReads.build_bowtie_index(correct_input), correct_output)
        with self.assertRaises(RuntimeError):
            shortReads.build_bowtie_index(empty_input)
    
    def test_align_reads(self):
        '''
        Test align_reads(reads_in)
        '''
        correct_input = 'Data/test_reads/reads_1.fq'
        correct_output = 'Data/test_reads/reads_1.sam'
        empty_input = ''
        self.assertEqual(shortReads.align_reads(correct_input), correct_output)
        with self.assertRaises(RuntimeError):
            shortReads.align_reads(empty_input)

if __name__ == '__main__':
    unittest.main()