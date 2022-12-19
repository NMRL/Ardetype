import unittest, sys, os
from pathlib import Path
sys.path.insert(0, os.path.dirname(os.path.dirname(Path(__file__).absolute())))
import subscripts.ardetype_utilities as au
import subscripts.src.modules as modules
from subscripts.ardetype_modules import Ardetype_module as am

class test_fun_1(unittest.TestCase):
    '''Testing the function using unittest library'''
    def test_valid_input(self):
        self.assertEqual(1, 1)

if __name__ == "__main__":
    unittest.main()