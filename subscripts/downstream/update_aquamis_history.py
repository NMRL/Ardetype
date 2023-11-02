import sys
sys.path.insert(0,'/mnt/beegfs2/home/groups/nmrl/bact_analysis/Ardetype/')
from subscripts.downstream import update_utilities as uu

#################
#Global variables
#################

full_path = uu.get_folder_path(__file__) #path to scripts
report_time = uu.get_current_timestamp()
arg_dict = {
    "aquamis"    : ["--aq", 'Full path to aquamis report (e.g. aquamis_report.csv)'],
    "batch"      : ["--bt", 'Batch name for the corresponding aquamis report (e.g. 2022-01-10-nmrl-qfxdna)']
}
proc_dict = uu.proc_dict

#defining custom report table processing parameters by changing specific default values
proc_dict["aq"] = proc_dict['default'].copy() #copy required to avoid overwriting the default
proc_dict["aq"]["delimiter"] = '\t'
proc_dict["aq"]["id_column"] = 'Sample_Name' #sample identifier
proc_dict["aq"]["primary_key"] = ['Sample_Name', 'analysis_batch_id']

#list of columns to ignore
ignore_list = ['batch']


##############
#Runtime logic
##############

if __name__ == '__main__':
    parser = uu.parse_arguments(arg_dict)
    proc_dict["aq"]["insert_batch"] = parser.bt
    uu.create_backup(full_path)
    current_tables = uu.find_current_tables(full_path, arg_dict, ignore_list=ignore_list)
    uu.update_files(arg_dict, parser, current_tables, report_time, full_path, proc_dict, ignore_list=ignore_list)
