import pandas as pd, time, os, sys, argparse
from pathlib import Path
#############################################################
#Placeholders for globals to be redefined in specific scripts
#############################################################
full_path = ''
report_time = ''
arg_dict = {}

##########
#Functions
##########
def get_folder_path(file_path:str):
    '''Returns full path to the folder where the file is located'''
    return str(Path(os.path.dirname(os.path.realpath(file_path))))

def get_current_timestamp(format:str="%m_%d_%Y"):
    '''Returns timestamp as string given timestamp format'''
    return time.strftime("%m_%d_%Y")

def parse_arguments(arg_dict) -> argparse.ArgumentParser:
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


def create_backup(full_path) -> None:
    '''
    Create a backup folder under global <full_path> if it does not exist
    '''
    if not os.path.isdir(f'{full_path}/backup/tables'):
        os.system(f'mkdir -p {full_path}/backup/tables')


def process_report(
    path_to_df:str, 
    delimiter:str=',', 
    col_rename_dict:dict={}, 
    batch_regex:list=[r'_[0-9]{8}_[0-9]{6}'],
    id_column:str='sample_id',
    id_regex:list=[r"_S[0-9]*"]
    ) -> pd.DataFrame:
    """
    Returns cleaned dataframe.
    path_to_df - path to file to read using pd.read_csv
    delimiter - delimiter to use (e.g. <,> for csv and <\t> for tsv)
    col_rename_dict - dictionary mapping initial and required column names
    batch_regex - list of strings representing regular expressions and substrings to replace in analysis_batch_id column
    id_column - name of the column that is expected to contain sample id
    id_regex - list of strings representing regular expressions and substrings to replace in id_columns
    """
    try:
        #Try to read dataframe
        df = pd.read_csv(path_to_df, sep=delimiter)

        #Rename columns if needed
        df = df.rename(columns=col_rename_dict)
        
        #To apply string operations with no surprises
        df[id_column] = df[id_column].astype(str)
        df.analysis_batch_id = df.analysis_batch_id.astype(str)
        #Remove substrings from id_column
        for regex in id_regex:
            df[id_column].str.replace(regex, '')

        #Remove substings from batch id
        for regex in batch_regex:
            df.analysis_batch_id = df.analysis_batch_id.str.replace(regex, '')

        return df
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def find_current_tables(full_path, arg_dict):
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


def update_history_file(
        path_to_new_table:str, 
        path_to_current_file:str, 
        current_combined_df:pd.DataFrame, 
        new_df:pd.DataFrame, 
        summary_name:str,
        report_time,
        full_path
        ):
    '''
    Given dataframe containing formatted batch of specific reports, 
    adds it to the combined dataframe, 
    moves previously generated combined table file to backup/, 
    moves new specific report tables to backyp/weekly/
    '''
    current_df = pd.concat([current_combined_df,new_df], sort=False)
    current_df.drop_duplicates(keep='last', inplace=True)
    os.system(f'mv {path_to_current_file} {full_path}/backup/')
    os.system(f'cp "{path_to_new_table}" {full_path}/backup/tables/{os.path.basename(path_to_new_table)}')
    current_df.to_csv(f'{full_path}/{summary_name}_{report_time}.csv', header=True, index=False)
    os.system(f'chmod 775 {full_path}/{summary_name}_{report_time}.csv')


def update_files(arg_dict:dict, parser, current_tables:dict, report_time:str, full_path:str):
    '''
    Wrapper function to update files in a loop.
    '''
    for report, arg in zip(arg_dict.keys(), vars(parser)):
            path_to_new_table    = getattr(parser,arg)
            path_to_current_file = current_tables[f'{report}_summary'][0]
            current_combined_df  = current_tables[f'{report}_summary'][1]
            new_df               = process_report(path_to_new_table)
            summary_name         = f'{report}_summary'
            update_history_file(
                path_to_new_table ,
                path_to_current_file,
                current_combined_df,
                new_df,
                summary_name,
                report_time,
                full_path
            )