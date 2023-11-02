import sys, os, time
from pathlib import Path
sys.path.insert(0,'/mnt/beegfs2/home/groups/nmrl/bact_analysis/Ardetype/')
from subscripts.downstream import update_utilities as uu

#################
#Global variables
#################

full_path = uu.get_folder_path(__file__) #path to scripts
report_time = uu.get_current_timestamp() #timestamp
arg_dict = {
    "kleborate"    : ["--kbt", "kleborate_report.csv"],
    "ectyper"      : ["--ect", "ectyper_report.csv"],
    "stecfinder"   : ["--stf", "stectfinder_report.csv"],
    "agrvate"      : ["--agr", "agrvate_report.csv"],
    "seqsero2"     : ["--ss2", "seqsero2_report.csv"],
    "sistr"        : ["--sst", "sistr_report.csv"],
    "lissero"      : ["--lss", "lissero_report.csv"],
    "meningotype"  : ["--mnt", "meningotype_report.csv"],
    "legsta"       : ["--lgt", "legsta_report.csv"],
    "chewbbaca_qc" : ["--cbc", "chewbbaca_qc_report.csv"],
    "lrefinder"    : ["--lrf", "lrefinder_report.csv"],
}

##############
#Runtime logic
##############

if __name__ == '__main__':
    parser = uu.parse_arguments(arg_dict)
    uu.create_backup(full_path)
    current_tables = uu.find_current_tables(full_path, arg_dict)
    uu.update_files(arg_dict, parser, current_tables, report_time, full_path)
