import unittest, pandas as pd, os
from subscripts.ardetype_modules import Ardetype_module as am


class test_module(unittest.TestCase):
    '''Testing methods of the module class'''


    #############################################################
    
    # Pre-testing configurations
    
    #############################################################


    # Methods used to verify dataframe equality based on https://stackoverflow.com/questions/38839402/how-to-use-assert-frame-equal-in-unittest
    def assertDataframeEqual(self, a, b, msg):
        try:
            pd.testing.assert_frame_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e


    def setUp(self):
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)


    @staticmethod
    def create_nested_dir_struct(branch_count:int=3, leave_count:int=3, file_count:int=1, root_name:str='top') -> list:
        '''
        Method is used to create nested folders with files to be used in testing process.
        The folder tree will have height of 3, which cannot be changed.
        root_name defines the name of the top-most folder
        branch_count - number of non-leave folders (containing subfolders)
        leave_count - number of leave folders (with no subfolders)
        file_count - number of files stored at each level
        
        Returns list of all created files.
        '''
        file_list = []
        for i in range(branch_count):
            for j in range(leave_count):
                os.makedirs(f'./{root_name}/middle{i+1}/bottom{j+1}', exist_ok=True)
                for _ in range(file_count):
                    open(f'./{root_name}/middle{i+1}/bottom{j+1}/{str(uuid.uuid4())}','a').close()
                    file_list.append(f'./{root_name}/middle{i+1}/bottom{j+1}/{str(uuid.uuid4())}')
        return file_list


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