import sys
sys.path.insert(0,'/home/group/pipelines/Ardetype/')
from subscripts.downstream import update_utilities as uu

#################
#Global variables
#################

full_path = uu.get_folder_path(__file__) #path to scripts
report_time = uu.get_current_timestamp()
arg_dict = {
    "k2contigs"  : ["--k2c", "Full path to the kraken2 contigs report (e.g. kraken2contigs_report.csv)"],
    "quast"      : ["--qst", "Full path to the quast report (e.g. quast_report.csv)"],
    "aquamis_qc" : ["--aqc", "Full path to the quast report (e.g. pointfinder_report.csv)"],
    "virfinder"  : ["--vir", "Full path to the quast report (e.g. virulencefinder_summary.csv)"]
}
proc_dict = uu.proc_dict
#setting primary key to None to ensure that records will be deduplicated only if all column value match
proc_dict['qst'] = proc_dict['default'].copy()
proc_dict['qst']['id_regex'].append('_quast')
proc_dict['k2c'] = proc_dict['default'].copy()
proc_dict['k2c']['id_regex'].append('_kraken2_contigs_report.txt')
proc_dict['vir'] = proc_dict['default'].copy()
proc_dict['k2c']['primary_key'] = None
proc_dict['vir']['primary_key'] = None

##############
#Runtime logic
##############

if __name__ == '__main__':
    parser = uu.parse_arguments(arg_dict)
    uu.create_backup(full_path)
    current_tables = uu.find_current_tables(full_path, arg_dict)
    uu.update_files(arg_dict, parser, current_tables, report_time, full_path, proc_dict)
