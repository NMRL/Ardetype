import sys
sys.path.insert(0,'/home/group/pipelines/Ardetype/')
from subscripts.downstream import update_utilities as uu

#################
#Global variables
#################

full_path = uu.get_folder_path(__file__) #path to scripts
report_time = uu.get_current_timestamp()
arg_dict = {
    "resistance"    : ["--hamr", 'Full path to the hamronization summary report (e.g. summarized_resistance_profile_(sequencing_date)_(batch).tsv)'],
    "resfinder"     : ["--rf", 'Full path to the aggregated resfinder phenotype report (e.g. resfinder_pheno_table_gathered.csv)'],
    "pointfinder"   : ["--pf", 'Full path to the aggregated pointfinder report (e.g. pointfinder_report.csv)'],
    "amrfinder_mut" : ["--afm", 'Full path to the aggregated amrfinder mutation report (e.g. amrfp_mutation_report.csv)'],

}
proc_dict = uu.proc_dict

#defining custom report table processing parameters by changing specific default values
proc_dict['hamr'] = proc_dict['default'].copy() #copy required to avoid overwriting the default
proc_dict['hamr']["delimiter"] = '\t' #to read table correctly
proc_dict['hamr']["id_column"] = 'input_file_name' #sample identifier
proc_dict['hamr']["id_regex"].append('.amr.alignment') #add more substrings to remove in id_column
proc_dict['hamr']["primary_key"] = None #if set to None, all columns will be used to infer uniqueness of rows
proc_dict['hamr']["history_ext"] = '.csv'

# resfinder custom preprocessing
proc_dict['rf'] = proc_dict['default'].copy() #copy required to avoid overwriting the default
proc_dict['rf']["id_regex"].append('_resfinder_pheno.txt') #add more substrings to remove in id_column

# amf+ custom preprocessing
proc_dict['afm'] = proc_dict['default'].copy() #copy required to avoid overwriting the default
proc_dict['afm']["id_regex"].append('_amrfinderplus_point.tab') #add more substrings to remove in id_column

#setting primary key to None to ensure that records will be deduplicated only if all column value match
for report in ['pf']:
    proc_dict[report] = proc_dict['default'].copy()
    proc_dict[report]['primary_key'] = None


##############
#Runtime logic
##############

if __name__ == '__main__':
    parser = uu.parse_arguments(arg_dict)
    uu.create_backup(full_path)
    current_tables = uu.find_current_tables(full_path, arg_dict)
    uu.update_files(arg_dict, parser, current_tables, report_time, full_path, proc_dict)
