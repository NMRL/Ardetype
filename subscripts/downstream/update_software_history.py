import sys
sys.path.insert(0,'/home/group/pipelines/Ardetype/')
from subscripts.downstream import update_utilities as uu

#################
#Global variables
#################

full_path = uu.get_folder_path(__file__) #path to scripts
report_time = uu.get_current_timestamp() #timestamp
#defining cmd args
arg_dict = {
    "software"   : ["--log",   "software_log.csv"]
}


proc_dict = uu.proc_dict
#defining custom report table processing parameters by changing specific default values
proc_dict["log"] = proc_dict['default'].copy() #copy required to avoid overwriting the default
proc_dict["log"]["id_column"] = 'analysis_batch_id' #sample identifier
proc_dict["log"]["primary_key"] = ["analysis_batch_id","tool","version"]

#list of columns to ignore
ignore_list = ['batch']


##############
#Runtime logic
##############

if __name__ == '__main__':
    parser = uu.parse_arguments(arg_dict)
    uu.create_backup(full_path)
    current_tables = uu.find_current_tables(full_path, arg_dict)
    uu.update_files(arg_dict, parser, current_tables, report_time, full_path, proc_dict)
