import unittest

class test_fun_1(unittest.TestCase):
    '''Testing the function using unittest library'''
    def test_valid_input(self):
        self.assertEqual(1, 1)

if __name__ == "__main__":
    unittest.main()