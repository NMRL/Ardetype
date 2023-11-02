import sys
sys.path.insert(0,'/mnt/beegfs2/home/groups/nmrl/bact_analysis/Ardetype/')
from subscripts.downstream import update_utilities as uu

#################
#Global variables
#################

full_path = uu.get_folder_path(__file__) #path to scripts
report_time = uu.get_current_timestamp()
arg_dict = {
    "plasmidfinder" : ["--plf", "Full path to the plasmidfinder summary report (e.g. plasmidfinder_summary.csv)"],
    "mobtyper"      : ["--mbt", "Full path to the aggregated mob_typer summary report (e.g. mobtyper_summary.csv)"]
}
proc_dict = uu.proc_dict
#setting primary key to None to ensure that records will be deduplicated only if all column value match
for report in ['plf', 'mbt']:
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
