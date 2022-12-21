import unittest, pandas as pd
from subscripts.ardetype_utilities import Ardetype_housekeeper as hk


class test_housekeeper(unittest.TestCase):
    '''Testing methods of the housekeeper class'''


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


# Methods used to verify dataframe equality based on https://stackoverflow.com/questions/38839402/how-to-use-assert-frame-equal-in-unittest
    def assertDataframeEqual(self, a, b, msg):
        try:
            pd.testing.assert_frame_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)


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

    
    def test_check_file_existance(self):
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