import os
import sys
import unittest
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from ectyper import (ectyper, loggingFunctions)

class TestPrediction(unittest.TestCase):
    
    def test_no_prediction(self):
        pass
    def test_unpaired_prediction(self):
        pass
    def test_one_prediction(self):
        pass
    def test_two_prediction(self):
        pass

if __name__ == '__main__':
    loggingFunctions.initialize_logging()
    unittest.main()
