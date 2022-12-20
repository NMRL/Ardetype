import unittest
import os
from subscripts.ardetype_utilities import Ardetype_housekeeper as hk
from subscripts.ardetype_modules import Ardetype_module as am

class test_housekeeper(unittest.TestCase):
    '''Testing methods of the housekeeper class'''
    def test_asign_perm_rec(self):
        self.assertEqual(1, 1)

'''
asign_perm_rec
check_file_existance
check_file_multiplicity
create_sample_sheet
edit_nested_dict
extract_log_id
filter_contigs_length
find_in_nested_dict
find_job_logs
get_all_keys
install_snakemake
map_new_column
name_job_logs
parse_arguments
parse_folder
parse_snakemake_log
printProgressBar
read_json_dict
read_yaml
remove_old_files
rename_file
type_contigs_api
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


'''
add_fasta_samples
add_module_targets
add_output_dir
add_taxonomy_column
check_job_completion
check_module_output
clear_working_directory
files_to_wd
fill_input_dict
fill_sample_sheet
fill_target_list
fold_output
make_output_dir
receive_sample_sheet
remove_invalid_samples
run_module
run_module_cluster
save_removed
set_permissions
submit_module_job
supply_sample_sheet
unfold_output
write_module_config
write_sample_sheet
'''


if __name__ == "__main__":
    unittest.main()