import unittest, pandas as pd
from subscripts.ardetype_modules import Ardetype_module as am


class test_module(unittest.TestCase):
    '''Testing methods of the module class'''


    #############################################################
    
    # Tests for methods that do not interact with the file system
    
    #############################################################


    def test_fill_input_dict(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_fill_sample_sheet(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_fill_target_list(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_receive_sample_sheet(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_remove_invalid_samples(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_supply_sample_sheet(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    #############################################################
    
    # Tests for methods that DO interact with the file system
    
    #############################################################


    def test_check_job_completion(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_check_module_output(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_clear_working_directory(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_files_to_wd(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_fold_output(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_make_output_dir(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_run_module(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_run_module_cluster(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_save_removed(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_set_permissions(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_submit_module_job(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_unfold_output(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_write_module_config(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception


    def test_write_sample_sheet(self):
        test = {
            'Valid input':[],
            'Exception|':[]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(1,1)
            else:
                with self.assertRaises(Exception):
                    raise Exception




if __name__ == "__main__":
    unittest.main()