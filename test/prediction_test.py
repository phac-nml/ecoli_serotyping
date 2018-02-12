from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
import datetime
import os
import sys
import unittest

import pandas as pd
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from ectyper import commandLineOptions, ectyper, predictionFunctions, definitions


class TestPrediction(unittest.TestCase):

    def test_no_prediction(self):
        print("\nTesting %s" % self.id())
        genome_file = 'test/Data/prediction/no_hit.fna'
        df, raw_df = self.helper(genome_file)
        self.assertEqual(df.shape[0], 1)
        self.assertEqual(df.O_prediction[0], '-')
        self.assertEqual(df.H_prediction[0], '-')
        self.assertTrue(raw_df.empty)

    def test_unpaired_prediction(self):
        print("\nTesting %s" % self.id())
        genome_file = 'test/Data/prediction/unpaired.fna'
        df, raw_df = self.helper(genome_file)
        self.assertEqual(df.shape[0], 1)
        self.assertEqual(df.O_prediction[0], '-')
        self.assertEqual(df.O_info[0], 'Only unpaired alignments found')
        self.assertEqual(df.H_prediction[0], '-')
        self.assertFalse(raw_df.empty)

    def test_lone_unpaired_prediction(self):
        print("\nTesting %s" % self.id())
        genome_file = 'test/Data/prediction/lone_unpaired.fna'
        df, raw_df = self.helper(genome_file)
        self.assertEqual(df.shape[0], 1)
        self.assertEqual(df.O_prediction[0], 'O64')
        self.assertEqual(df.O_info[0], 'Lone unpaired alignment found')
        self.assertEqual(df.H_prediction[0], '-')
        self.assertFalse(raw_df.empty)

    def test_one_prediction(self):
        genome_file = 'test/Data/prediction/one_hit.fna'
        df, raw_df = self.helper(genome_file)
        self.assertEqual(df.shape[0], 1)
        self.assertEqual(df.O_prediction[0], '-')
        self.assertEqual(df.H_prediction[0], 'H7')
        self.assertFalse(raw_df.empty)

    def test_two_prediction(self):
        print("\nTesting %s" % self.id())
        genome_file = 'test/Data/prediction/two_hit.fna'
        df, raw_df = self.helper(genome_file)
        self.assertEqual(df.shape[0], 1)
        self.assertEqual(df.O_prediction[0], 'O157')
        self.assertEqual(df.H_prediction[0], 'H1')
        self.assertFalse(raw_df.empty)

    def helper(self, genome_file):
        args = commandLineOptions.parse_command_line(['-i', 'test'])
        genome_file = genome_file
        with definitions.TEMPDIR() as temp_dir:
            prediction_file = os.path.join(temp_dir + 'output.csv')
            ectyper.run_prediction([genome_file], args, prediction_file)
            basename, extension = os.path.splitext(prediction_file)
            raw_prediction_file = ''.join([basename, '_raw', extension])
            predictionFunctions.add_non_predicted(
                [os.path.splitext(os.path.basename(genome_file))[0]], prediction_file)
            return pd.read_csv(prediction_file, encoding='utf-8'), pd.read_csv(raw_prediction_file, encoding='utf-8')


if __name__ == '__main__':
    unittest.main()
