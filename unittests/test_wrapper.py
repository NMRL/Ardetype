import unittest
from subscripts.ardetype_utilities import Ardetype_housekeeper as hk
from subscripts.ardetype_modules import Ardetype_module as am



class test_housekeeper(unittest.TestCase):
    '''Testing the function using unittest library'''
    def test_valid_input(self):
        self.assertEqual(1, 1)

if __name__ == "__main__":
    unittest.main()