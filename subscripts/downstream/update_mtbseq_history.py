import pandas as pd, time, os, sys, argparse
from pathlib import Path
pd.set_option('mode.chained_assignment', None)

#################
#Global variables
#################

full_path = str(Path(os.path.dirname(os.path.realpath(__file__)))) #path to scripts
report_time = time.strftime("%m_%d_%Y") #timestamp
arg_dict = {
    "variant"         : ["--var" , "MTBseq_variant_report.csv"],
    "statistics"      : ["--stat", "MTBseq_statistics_report.csv"],
    "classification"  : ["--cls" , "MTBseq_classification_report.csv"],
}

##########
#Functions
##########

def parse_arguments() -> argparse.ArgumentParser:
    '''Parsing predefined command-line arguments and return argparse class object.'''

    parser = argparse.ArgumentParser(description='A script to update database file with new batch data.') 
    
    for arg in arg_dict.keys():
        parser.add_argument(
            arg_dict[arg][0],
            f'--{arg}',
            metavar='\b',
            help = f'Full path to the {arg} report (e.g. {arg_dict[arg][1]})',
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


def process_report(path_to_df:str, delimiter:str='.csv') -> pd.DataFrame:
    """ 
    Returns cleaned dataframe. 
    """

    if path_to_df.endswith(delimiter):
        try:
            df = pd.read_csv(path_to_df)
            df['analysis_batch_id'] = df['analysis_batch_id'].str.replace(r'_[0-9]{8}_[0-9]{6}', '')
            return df
        except pd.errors.EmptyDataError:
            return pd.DataFrame()
    else:
        raise ValueError(f'{path_to_df} is expected to end with {delimiter}')


def find_current_tables():
    '''
    Checks files in the folder where the script is placed. 
    Returns tuple containing pd.Dataframe object generated from that file and path to the file, if file name contains "resistance_summary" substring. 
    Returns None if file is not found.
    '''
    current_dict = {}
    for file in os.listdir(full_path):
        path_to_file = f'{full_path}/{file}'
        for arg in arg_dict:
            if os.path.isfile(path_to_file) and f'{arg}_summary' in file:
                try:
                    current_dict[f'{arg}_summary'] = [
                        path_to_file, 
                        pd.read_csv(path_to_file, low_memory=False)
                    ]
                except pd.errors.EmptyDataError:
                    current_dict[f'{arg}_summary'] = [
                        path_to_file, 
                        pd.DataFrame()
                    ]
    return current_dict


def update_combined_agnostic_file(
        path_to_new_table:str, 
        path_to_current_file:str, 
        current_combined_df:pd.DataFrame, 
        new_df:pd.DataFrame, 
        summary_name:str
        ):
    '''
    Given dataframe containing formatted batch of specific reports, 
    adds it to the combined dataframe, 
    moves previously generated combined table file to backup/, 
    moves new specific report tables to backyp/weekly/
    '''
    current_df = pd.concat([current_combined_df,new_df], sort=False)
    current_df.drop_duplicates(keep='last', inplace=True)
    # print(current_df, len(current_df), path_to_current_file)
    os.system(f'mv {path_to_current_file} {full_path}/backup/')
    os.system(f'cp "{path_to_new_table}" {full_path}/backup/tables/{os.path.basename(path_to_new_table)}')
    current_df.to_csv(f'{full_path}/{summary_name}_{report_time}.csv', header=True, index=False)
    os.system(f'chmod 775 {full_path}/{summary_name}_{report_time}.csv')


##############
#Runtime logic
##############

if __name__ == '__main__':
    #copy current summary table to backup
    create_backup()
    parser                 = parse_arguments()
    #read current summary table
    current_tables = find_current_tables() 
    # print(len(current_tables), len(vars(parser)), len(arg_dict.keys()))
    #read and preprocess new tables
    for report, arg in zip(arg_dict.keys(), vars(parser)):
        path_to_new_table    = getattr(parser,arg)
        path_to_current_file = current_tables[f'{report}_summary'][0]
        current_combined_df  = current_tables[f'{report}_summary'][1]
        new_df               = process_report(path_to_new_table)
        summary_name         = f'{report}_summary'
        # print(path_to_new_table, path_to_current_file, summary_name)
        update_combined_agnostic_file(
            path_to_new_table ,
            path_to_current_file,
            current_combined_df,
            new_df,
            summary_name
        )

