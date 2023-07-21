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
        '-p', 
        '--pipeline', 
        metavar='\b', 
        help = 'Full path to aquamis report (e.g. aquamis_report.csv)', 
        default=None, 
        required=True
        )

    parser.add_argument(
        '-b', 
        '--batch', 
        metavar='\b', 
        help = 'Batch name for the corresponding aquamis report (e.g. 2022-01-10-nmrl-qfxdna)', 
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


def process_report(report_path:str, parser) -> pd.DataFrame:
    '''Reads pipe_report_path and db_file_path into pandas dataframes. Returns tuple of two dataframes.'''
    df = pd.read_csv(report_path, sep='\t') #READ PIPELINE REPORT
    df.insert(1,'analysis_batch_id',[str(parser.batch) for _ in df.index])
    return df


def find_current_table():
    '''
    Checks files in the folder where the script is placed. 
    Returns tuple containing pd.Dataframe object generated from that file and path to the file, if file name contains "aquamis_summary" substring. 
    Returns None if file is not found.
    '''
    aquamis_summary, aqm_df      = '', None
    for file in os.listdir(full_path):
        path_to_file = f'{full_path}/{file}'
        if os.path.isfile(path_to_file) and 'aquamis_summary' in file: 
            aquamis_summary, aqm_df      = path_to_file, pd.read_csv(path_to_file, low_memory=False)
    if aqm_df is not None:
        return aquamis_summary, aqm_df
    else:
        return None


def update_aquamis_file(path_to_new_table:str, path_to_current_file:str, current_combined_df:pd.DataFrame, new_df:pd.DataFrame):
    '''
    Given dataframe containing formatted batch of aquamis, 
    adds it to the combined dataframe, 
    moves previously generated combined table file to backup/, 
    moves new resistance table to backyp/weekly/,
    saves updated combined dataframe as "aquamis_summary" file with current timestamp and sets its permissions to 775.
    '''
    current_df = pd.concat([current_combined_df,new_df], sort=False)
    current_df = current_df.astype(str)
    current_df.drop_duplicates(subset=['Sample_Name', 'analysis_batch_id'], keep='last', inplace=True)
    os.system(f'mv {path_to_current_file} {full_path}/backup/')
    os.system(f'cp "{path_to_new_table}" {full_path}/backup/tables/{os.path.basename(path_to_new_table)}')
    current_df.to_csv(f'{full_path}/aquamis_summary_{report_time}.csv', header=True, index=False)
    os.system(f'chmod 775 {full_path}/aquamis_summary_{report_time}.csv')


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
    aqm_df = process_report(os.path.abspath(parser.pipeline), parser) 
    #update hamronization-aggregated summary table
    update_aquamis_file( 
        path_to_new_table=parser.pipeline ,
        current_combined_df=current_combined_tuple[1], 
        path_to_current_file=current_combined_tuple[0],
        new_df=aqm_df
    )
