import sys
sys.path.insert(0,'/home/group/pipelines/Ardetype/')
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
    "spatyper"     : ["--spa", "spatyper_report.csv"],
    "shigatyper"   : ["--sht", "shigatyper_report.csv"],
    "seroba"       : ["--sba", "seroba_report.csv"],
    "emmtyper"     : ["--emm", "emmtyper_report.csv"],
    "hicap"        : ["--hic", "hicap_report.csv"],
}
proc_dict = uu.proc_dict

#setting primary key to None to ensure that records will be deduplicated only if all column value match
for report in ['lrf', 'spa']:
    proc_dict[report] = proc_dict['default'].copy()
    proc_dict[report]['primary_key'] = None #infer uniqueness of rows from all columns

#setting custom id_column for kleborate report
proc_dict['kbt'] = proc_dict['default'].copy()
proc_dict['kbt']['col_rename_dict'] = {'strain':'sample_id'}

#custom id_regex for chewbacca qc 
proc_dict['cbc'] = proc_dict['default'].copy()
proc_dict['cbc']['id_regex'].append('_contigs')

#custom id_regex for stecfinder 
proc_dict['stf'] = proc_dict['default'].copy()
proc_dict['stf']['id_regex'].append('_bact_reads_classified')

#custom id_regex for lrefinder
proc_dict['lrf'] = proc_dict['default'].copy()
proc_dict['lrf']['id_regex'].append('.pos')


##############
#Runtime logic
##############

if __name__ == '__main__':
    parser = uu.parse_arguments(arg_dict)
    uu.create_backup(full_path)
    current_tables = uu.find_current_tables(full_path, arg_dict)
    uu.update_files(arg_dict, parser, current_tables, report_time, full_path, proc_dict)
