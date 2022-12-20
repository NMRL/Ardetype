import unittest
import os
from subscripts.ardetype_utilities import Ardetype_housekeeper as hk
from subscripts.ardetype_modules import Ardetype_module as am



class test_housekeeper(unittest.TestCase):
    '''Testing methods of the housekeeper class'''
    def test_edit_nested_dict(self):
        test = {
            'Valid case':[{'level_1':{'level_2':{'level_3':1}}}, 'level_3', 0, 0],
            'Missing parameter':[{'level_1':{'level_2':{'level_3':1}}}, 'level_4', 0, None],
            'Wrong nested iterable format':[[], '', '', AttributeError],
            'Wrong attribute format':[{}, ['XYZ'], '', TypeError]
        }
        for case in test:
            if 'format' not in case:
                self.assertEqual(hk.edit_nested_dict(test[case][0],test[case][1],test[case][2]),test[case][3])
            else:
                with self.assertRaises(test[case][3]):
                    hk.edit_nested_dict(test[case][0],test[case][1],test[case][2])

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

#Methods that operate on python objects
'''
edit_nested_dict
find_in_nested_dict
get_all_keys
map_new_column
'''


class test_module(unittest.TestCase):
    '''Testing methods of the module class'''
    def test_valid_input(self):
        self.assertEqual(1, 1)


#Methods that operate on file system
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
#Methods that operate on python objects

'''
fill_input_dict
fill_sample_sheet
fill_target_list
receive_sample_sheet
remove_invalid_samples
supply_sample_sheet
'''


if __name__ == "__main__":
    unittest.main()