import pandas as pd, time, os, sys, argparse
from pathlib import Path
pd.set_option('mode.chained_assignment', None)

#################
#Global variables
#################

full_path = str(Path(os.path.dirname(os.path.realpath(__file__)))) #path to scripts
report_time = time.strftime("%m_%d_%Y") #timestamp

##########
#Functions
##########

def parse_arguments() -> argparse.ArgumentParser:
    '''Parsing predefined command-line arguments and return argparse class object.'''

    parser = argparse.ArgumentParser(description='A script to update database file with new batch data.') 
    parser.add_argument(
        '-k2r', 
        '--kraken2_reads', 
        metavar='\b', 
        help = 'Full path to the kraken2 reads report (e.g. kraken2reads_report.csv)', 
        default=None, 
        required=True
        )
    parser.add_argument(
        '-k2c', 
        '--kraken2_contigs', 
        metavar='\b', 
        help = f'Full path to the kraken2 contigs report (e.g. kraken2contigs_report.csv)', 
        default=None, 
        required=True
        )
    parser.add_argument(
        '-qst', 
        '--quast', 
        metavar='\b', 
        help = f'Full path to the quast report (e.g. pointfinder_report.csv)', 
        default=None, 
        required=True
        )


    #if script is run without arguments - print help message & exit
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    return parser.parse_args()


def create_backup() -> None:
    '''
    Create a backup folder under  if it does not exist
    '''
    if not os.path.isdir(f'{full_path}/backup/tables'):
        os.system(f'mkdir -p {full_path}/backup/tables')


def process_report(path_to_df:str) -> pd.DataFrame:
    """ 
    Returns cleaned dataframe. 
    """

    if path_to_df.endswith(".csv"):
        df                 = pd.read_csv(path_to_df)
        df.sample_id = df.sample_id.astype(str)
        df.sample_id = df.sample_id.str.replace(r"_S[0-9]*", '')
        return df
    else:
        raise ValueError(f'{path_to_df} is expected to end with .csv')


def find_current_table():
    '''
    Checks files in the folder where the script is placed. 
    Returns tuple containing pd.Dataframe object generated from that file and path to the file, if file name contains "resistance_summary" substring. 
    Returns None if file is not found.
    '''
    k2reads_summary, k2r   = '', None
    k2contigs_summary, k2c = '', None
    for file in os.listdir(full_path):
        path_to_file = f'{full_path}/{file}'
        if os.path.isfile(path_to_file) and 'k2reads_summary' in file: 
            try:
                k2reads_summary, k2r = path_to_file, pd.read_csv(path_to_file, low_memory=False, sep='\t')
            except pd.errors.EmptyDataError:
                k2reads_summary, k2r = path_to_file, pd.DataFrame()
        elif os.path.isfile(path_to_file) and 'k2contigs_summary' in file:
            try:
                k2contigs_summary, k2c  = path_to_file, pd.read_csv(path_to_file, low_memory=False)
            except pd.errors.EmptyDataError:
                k2contigs_summary, k2c = path_to_file, pd.DataFrame()
        elif os.path.isfile(path_to_file) and 'quast_summary' in file:
            try:
                quast_summary, qst     = path_to_file, pd.read_csv(path_to_file, low_memory=False)
            except pd.errors.EmptyDataError:
                quast_summary, qst     = path_to_file, pd.DataFrame()
    if k2r is not None and k2c is not None and qst is not None:
        return k2reads_summary, k2r, k2contigs_summary, k2c, quast_summary, qst
    else:
        return None


def update_combined_agnostic_file(path_to_new_table:str, path_to_current_file:str, current_combined_df:pd.DataFrame, new_df:pd.DataFrame, type='kraken2_reads'):
    '''
    Given dataframe containing formatted batch of resistance, 
    adds it to the combined dataframe, 
    moves previously generated combined table file to backup/, 
    moves new resistance table to backyp/weekly/,
    saves updated combined dataframe as "resistance_summary" file with current timestamp and sets its permissions to 775.
    '''
    if type == 'kraken2_reads':
        current_df = pd.concat([current_combined_df,new_df], sort=False)
        current_df.drop_duplicates(keep='last', inplace=True)
        os.system(f'mv {path_to_current_file} {full_path}/backup/')
        os.system(f'cp "{path_to_new_table}" {full_path}/backup/tables/{os.path.basename(path_to_new_table)}')
        current_df.to_csv(f'{full_path}/k2reads_summary_{report_time}.csv', header=True, index=False)
        os.system(f'chmod 775 {full_path}/k2reads_summary_{report_time}.csv')
    elif type == 'kraken2_contigs':
        #separate contam.records
        current_df = pd.concat([current_combined_df,new_df], sort=False)
        current_df.drop_duplicates(keep='last', inplace=True)
        os.system(f'mv {path_to_current_file} {full_path}/backup/')
        os.system(f'cp "{path_to_new_table}" {full_path}/backup/tables/{os.path.basename(path_to_new_table)}')
        current_df.to_csv(f'{full_path}/k2contigs_summary_{report_time}.csv', header=True, index=False)
        os.system(f'chmod 775 {full_path}/k2contigs_summary_{report_time}.csv')
    elif type == 'quast':
        current_df = pd.concat([current_combined_df,new_df], sort=False)
        current_df.drop_duplicates(keep='last', inplace=True)
        os.system(f'mv {path_to_current_file} {full_path}/backup/')
        os.system(f'cp "{path_to_new_table}" {full_path}/backup/tables/{os.path.basename(path_to_new_table)}')
        current_df.to_csv(f'{full_path}/quast_summary_{report_time}.csv', header=True, index=False)
        os.system(f'chmod 775 {full_path}/quast_summary_{report_time}.csv')


##############
#Runtime logic
##############

if __name__ == '__main__':
    #copy current summary table to backup
    create_backup()
    parser                 = parse_arguments()
    #read current summary table
    current_combined_tuple = find_current_table() 
    #read and preprocess new tables
    k2r = process_report(os.path.abspath(parser.kraken2_reads)) 
    k2c = pd.read_csv(os.path.abspath(parser.kraken2_contigs))
    qst  = pd.read_csv(os.path.abspath(parser.quast))
    #update hamronization-aggregated summary table
    update_combined_agnostic_file( 
        path_to_new_table=parser.kraken2_reads,
        current_combined_df=current_combined_tuple[1], 
        path_to_current_file=current_combined_tuple[0],
        new_df=k2r
    )

    #update resfinder-phenotype summary table
    update_combined_agnostic_file(
        path_to_new_table=parser.kraken2_contigs, 
        path_to_current_file=current_combined_tuple[2], 
        current_combined_df=current_combined_tuple[3], 
        new_df=k2c,
        type="kraken2_contigs"
    )

    #update pointfinder summary table
    update_combined_agnostic_file(
        path_to_new_table=parser.quast, 
        path_to_current_file=current_combined_tuple[4], 
        current_combined_df=current_combined_tuple[5], 
        new_df=qst,
        type="quast"
    )
