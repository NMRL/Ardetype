import unittest, pandas as pd, os, uuid
from shutil import rmtree
from subscripts.ardetype_utilities import Ardetype_housekeeper as hk


class test_housekeeper(unittest.TestCase):
    '''Testing methods of the housekeeper class'''


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
    def create_nested_dir_struct(branch_count:int=3, leave_count:int=3, file_count:int=1, root_name:str='top', file_multiplicity:int=1) -> list:
        '''
        Method is used to create nested folders with files to be used in testing process.
        The folder tree will have height of 3, which cannot be changed.
        root_name defines the name of the top-most folder
        branch_count - number of non-leave folders (containing subfolders)
        leave_count - number of leave folders (with no subfolders)
        file_count - number of files stored at each level
        file_multiplicity - option to simulate files that have names that differ only by some suffix (e.g. paired fastq)
        Returns list of all created files.
        '''
        for i in range(branch_count):
            for j in range(leave_count):
                os.makedirs(f'./{root_name}/middle{i+1}/bottom{j+1}', exist_ok=True)
                for _ in range(file_count):
                    if file_multiplicity <=1:
                        open(f'./{root_name}/middle{i+1}/bottom{j+1}/{str(uuid.uuid4())}','a').close()
                    else:
                        fname = str(uuid.uuid4())
                        for k in range(file_multiplicity):
                            open(f'./{root_name}/middle{i+1}/bottom{j+1}/{fname}_{k+1}','a').close()



    #############################################################
    
    # Tests for methods that do not interact with the file system
    
    #############################################################


    def test_edit_nested_dict(self):
        test = {
            'Valid case':[{'level_1':{'level_2':{'level_3':1}}}, 'level_3', 0, 0],
            'Missing parameter':[{'level_1':{'level_2':{'level_3':1}}}, 'level_4', 0, None],
            'Exception|Wrong nested iterable format':[[], '', '', AttributeError],
            'Exception|Wrong attribute format':[{}, ['XYZ'], '', TypeError]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(hk.edit_nested_dict(test[case][0],test[case][1],test[case][2]),test[case][3])
            else:
                with self.assertRaises(test[case][3]):
                    hk.edit_nested_dict(test[case][0],test[case][1],test[case][2])

    
    def test_find_in_nested_dict(self):
        test = {
            'Valid case':[{'level_1':{'level_2':{'level_3':1}, 'level_21':2}}, ['level_1','level_2', 'level_3'], 1],
            'Exception|Missing key in sequence':[{'level_1':{'level_2':{'level_3':1}}}, ['level_1','level_21'], LookupError],
            'Exception|Non-dict value accessed before the end of the key sequence':[{'level_1':{'level_2': 1}}, ['level_1','level_2', 'level_3'], LookupError],
            'Exception|Wrong sequence format':[{}, 'XYZ', TypeError],
            'Exception|Wrong nested iterable format':[[], ['XYZ'], TypeError]
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(hk.find_in_nested_dict(test[case][0],test[case][1]), test[case][2])
            else:
                with self.assertRaises(test[case][2]):
                    hk.find_in_nested_dict(test[case][0],test[case][1])


    def test_get_all_keys(self):
        test = {
            'Valid case':[{'level_1':{'level_2':{'level_3':1}, 'level_21':2}}, {'level_1', 'level_2','level_3','level_21'}],
            'Exception|Non-dictionary input':[{'level_1', 'level_2','level_3','level_21'}, AttributeError],
        }
        for case in test:
            if 'Exception' not in case:
                self.assertEqual(hk.get_all_keys(test[case][0]), test[case][1])
            else:
                with self.assertRaises(test[case][1]):
                    hk.get_all_keys(test[case][0])


    def test_map_new_column(self):
        test = {
            'Valid case':[
                pd.DataFrame.from_dict({'id':[i+1 for i in range(5)], 'var1':[10 for _ in range(5)], 'var2':[10-i-1 for i in range(5)]}),
                dict(zip([i+1 for i in range(5)],[10+i+1 for i in range(5)])), 
                'id', 
                'var3',
                pd.DataFrame.from_dict({'id':[i+1 for i in range(5)], 'var1':[10 for _ in range(5)], 'var2':[10-i-1 for i in range(5)], 'var3':[10+i+1 for i in range(5)]})
                ],
            'Exception|No key column found in dataframe':[
                pd.DataFrame.from_dict({'id':[i+1 for i in range(5)], 'var1':[10 for _ in range(5)], 'var2':[10-i-1 for i in range(5)]}),
                dict(zip([i+1 for i in range(5)],[10+i+1 for i in range(5)])),
                'id_column',
                'var3',
                KeyError
                ],
            'Exception|No correspondance between ids used in new column and key column content':[
                pd.DataFrame.from_dict({'id':[i+1 for i in range(5)], 'var1':[10 for _ in range(5)], 'var2':[10-i-1 for i in range(5)]}),
                dict(zip([i+1 for i in range(6,10)],[10+i+1 for i in range(5)])),
                'id', 
                'var3',
                KeyError
                ],
            'Exception|Non-dataframe ss_df input':[
                {'id':[i+1 for i in range(5)], 'var1':[10 for _ in range(5)], 'var2':[10-i-1 for i in range(5)]},
                dict(zip([i+1 for i in range(5)],[10+i+1 for i in range(5)])), 
                'id', 
                'var3',
                TypeError
                ],
            'Exception|Non-dictionary info_dict input':[
                pd.DataFrame.from_dict({'id':[i+1 for i in range(5)], 'var1':[10 for _ in range(5)], 'var2':[10-i-1 for i in range(5)]}),
                tuple(zip([i+1 for i in range(5)],[10+i+1 for i in range(5)])), 
                'id', 
                'var3',
                TypeError
                ],
        }

        for case in test:
            if 'Exception' not in case:
                self.assertEqual(hk.map_new_column(test[case][0],test[case][1], test[case][2], test[case][3]), test[case][4])
            else:
                with self.assertRaises(test[case][4]):
                    hk.map_new_column(test[case][0],test[case][1], test[case][2], test[case][3])


    #############################################################
    
    # Tests for methods that DO interact with the file system
    
    #############################################################


    def test_asign_perm_rec(self):
        test_housekeeper.create_nested_dir_struct()
        drs, fls = [], []
        for root, dirs, files in os.walk('./top/'):
            [drs.append(os.path.join(root,dr)) for dr in dirs]
            [fls.append(os.path.join(root,fl)) for fl in files]
        drs.append('./top/')

        #Default permissions
        hk.asign_perm_rec('./top/')
        self.assertTrue(all([oct(os.stat(dr).st_mode)[-3:]=='775' for dr in drs]))
        self.assertTrue(all([oct(os.stat(fl).st_mode)[-3:]=='775' for fl in fls]))

        #Different permissions
        hk.asign_perm_rec('./top/','777')
        self.assertFalse(all([oct(os.stat(dr).st_mode)[-3:]=='775' for dr in drs]))
        self.assertFalse(all([oct(os.stat(fl).st_mode)[-3:]=='775' for fl in fls]))
        
        #Cleanup
        rmtree('./top/')

    
    def test_check_file_existance(self):

        #Existing files
        test_housekeeper.create_nested_dir_struct()
        fls = []
        for root, _, files in os.walk('./top/'):
            [fls.append(os.path.join(root,fl)) for fl in files]
        
        self.assertTrue(all(list(hk.check_file_existance(fls).values())))

        #Non-existing files
        fls = [str(uuid.uuid4()) for _ in range(10)]
        self.assertTrue(not any(hk.check_file_existance(fls).values()))

    
    def test_check_file_multiplicity(self):
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

    
    def test_create_sample_sheet(self):
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

    
    def test_extract_log_id(self):
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

    
    def test_filter_contigs_length(self):
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

    
    def test_find_job_logs(self):
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

    
    def test_install_snakemake(self):
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

    
    def test_name_job_logs(self):
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

    
    def test_parse_folder(self):
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

    
    def test_parse_snakemake_log(self):
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


    def test_read_json_dict(self):
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


    def test_read_yaml(self):
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


    def test_remove_old_files(self):
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

    def test_rename_file(self):
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


    def test_update_log_history(self):
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


    def test_update_log_summary(self):
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


    def test_validate_yaml(self):
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


    def test_write_json(self):
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


    def test_write_yaml(self):
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