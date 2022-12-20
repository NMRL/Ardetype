import unittest, pandas as pd
from subscripts.ardetype_utilities import Ardetype_housekeeper as hk
from subscripts.ardetype_modules import Ardetype_module as am


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


#Methods that operate on file system
'''
asign_perm_rec
check_file_existance
check_file_multiplicity
create_sample_sheet
extract_log_id
filter_contigs_length
find_job_logs
install_snakemake
name_job_logs
parse_folder
parse_snakemake_log
read_json_dict
read_yaml
remove_old_files
rename_file
update_log_history
update_log_summary
validate_yaml
write_json
write_yaml
'''




class test_module(unittest.TestCase):
    '''Testing methods of the module class'''
    def test_valid_input(self):
        self.assertEqual(1, 1)


    #############################################################
    
    # Tests for methods that do not interact with the file system
    
    #############################################################

'''
fill_input_dict
fill_sample_sheet
fill_target_list
receive_sample_sheet
remove_invalid_samples
supply_sample_sheet
'''
    #############################################################
    
    # Tests for methods that DO interact with the file system
    
    #############################################################

'''
check_job_completion
check_module_output
clear_working_directory
files_to_wd
fold_output
make_output_dir
run_module
run_module_cluster
save_removed
set_permissions
submit_module_job
unfold_output
write_module_config
write_sample_sheet
'''



if __name__ == "__main__":
    unittest.main()