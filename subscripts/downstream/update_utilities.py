import pandas as pd, time, os, sys, argparse
from pathlib import Path

#############################################################
#Placeholders for globals to be redefined in specific scripts
#############################################################

full_path = ''
report_time = ''
arg_dict = {}
proc_dict = {
    'default' : {
        "delimiter"       :',', 
        "col_rename_dict" : {}, 
        "batch_regex"     : [r'_[0-9]{8}_[0-9]{6}'],
        "id_column"       : 'sample_id',
        "id_regex"        : [r"_S[0-9]*"],
        "batch_column"    : "analysis_batch_id",
        'insert_batch'    : None,
        'primary_key'     : ['sample_id', "analysis_batch_id"],
        'history_ext'     : '.csv',
        'custom_functions': None
    }
}


##########
#Functions
##########
def get_folder_path(file_path:str):
    '''Returns full path to the folder where the file is located'''
    return str(Path(os.path.dirname(os.path.realpath(file_path))))

def get_current_timestamp(stamp_format:str="%m_%d_%Y"):
    '''Returns timestamp as string given timestamp format'''
    return time.strftime(stamp_format)

def parse_arguments(arg_dict:dict) -> argparse.ArgumentParser:
    '''Parsing predefined command-line arguments and return argparse class object.
    Expected input format:
    arg_dict = {
        "<report_prefix>"    : ["--<abbreviation>", "Example help message specifying the expected file name, e.g. <report_prefix>_report.csv"],
    }
    '''
    parser = argparse.ArgumentParser(description='A script to update database file with new batch data.') 
    for arg in arg_dict.keys():
        parser.add_argument(
            arg_dict[arg][0],
            f'--{arg}',
            metavar='\b',
            help = arg_dict[arg][1],
            default=None, 
            required=True
        )

    #if script is run without arguments - print help message & exit
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    return parser.parse_args()


def create_backup(full_path:str) -> None:
    '''
    Creates a backup folder under global <full_path> if it does not exist
    '''
    if not os.path.isdir(f'{full_path}/backup/tables'):
        os.system(f'mkdir -p {full_path}/backup/tables')


def process_report(
        path_to_df:str, 
        delimiter:str=',', 
        col_rename_dict:dict={}, 
        batch_regex:list=[r'_[0-9]{8}_[0-9]{6}'],
        id_column:str='sample_id',
        id_regex:list=[r"_S[0-9]*"],
        insert_batch:str=None,
        batch_column:str='analysis_batch_id', 
        custom_functions:list=None
    ) -> pd.DataFrame:
    """
    path_to_df - path to file to read using pd.read_csv;
    delimiter - delimiter to use (e.g. <,> for csv and <\t> for tsv);
    col_rename_dict - dictionary mapping initial and required column names;
    batch_regex - list of strings representing regular expressions and substrings to replace in analysis_batch_id column;
    id_column - name of the column that is expected to contain sample id;
    id_regex - list of strings representing regular expressions and substrings to replace in id_columns;
    insert_batch - in order to allow supplying batch id through command line (see update_aquamis_history.py);
    batch_column - name of the batch column in case it is different from analysis_batch_id for some reason;
    custom_functions - list containing references to callable python functions to perform case-specific processing (see update_ardetype_history.py);
    """
    try:
        #Try to read dataframe
        df = pd.read_csv(path_to_df, sep=delimiter)

        #Rename columns if needed
        df = df.rename(columns=col_rename_dict)

        #Insert analysis_batch_id column with value supplied via command-line argument
        if insert_batch is not None:
            df.insert(1,batch_column,[insert_batch for _ in df.index])
        
        #To apply string operations with no surprises
        df[id_column] = df[id_column].astype(str)
        df[batch_column] = df[batch_column].astype(str)
        #Remove substrings from id_column
        for regex in id_regex:
            df[id_column].str.replace(regex, '')

        #Remove substings from batch id
        for regex in batch_regex:
            df.analysis_batch_id = df.analysis_batch_id.str.replace(regex, '')

        #Applying custom processing functions
        for func in custom_functions:
            df = func(df)
        
        return df
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def find_current_tables(full_path:str, arg_dict:dict, ignore_list:list=[]):
    '''
    Checks files in the folder where the script is placed. 
    Returns tuple containing pd.Dataframe object generated from that file and path to the file, if file name contains "<report_prefix>_summary" substring. 
    arg_dict - refers to arg_dict global variable redefind in specific module;
    full_path - refers to full_path global variable redefind in specific module;
    ignore_list - list of command-line arguments that do not refer to paths to file (e.g. batch number in update_aquamis_history.py);
    '''
    current_dict = {}
    for file in os.listdir(full_path):
        path_to_file = f'{full_path}/{file}'
        for arg in arg_dict:
            if arg not in ignore_list:
                if os.path.isfile(path_to_file) and f'{arg}_summary' in file:
                    try:
                        if file.endswith('.csv'):
                            current_dict[f'{arg}_summary'] = [
                                path_to_file, 
                                pd.read_csv(path_to_file, low_memory=False)
                            ]
                        elif file.endswith('.tsv'):
                            current_dict[f'{arg}_summary'] = [
                                path_to_file, 
                                pd.read_csv(path_to_file, low_memory=False, delimiter='\t')
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
        report_time:str,
        full_path:str,
        primary_key:list,
        history_ext
        ):
    '''
    Given dataframe containing formatted batch of specific reports,
    adds it to the combined dataframe;
    removes duplicates either by using primary_key or entire row as unique row id;
    moves previously generated combined table file to <full_path>backup/;
    moves new specific report tables to <full_path>backup/weekly/;
    applies 775 permission to the new file.

    path_to_new_table - full path to the new report file (to be backed-up);
    path_to_current_file - full path to the current history file (to be backed-up);
    current_combined_df - dataframe corresponding to the current history file;
    new_df - dataframe corresponding to the new report file;
    summary_name - summary file prefix expected as <report_prefix>_summary string;
    report_time - refers to report_time global variable redefind in specific module (timestamp string expected);
    full_path - refers to full_path global variable redefind in specific module;
    primary_key - None or list of strings representing combination of column names to be used as unique row identifier for the history file;
    history_ext - string representing extension of the history file (e.g. `.csv`);
    '''

    current_df = pd.concat([current_combined_df,new_df], sort=False)
    if primary_key is not None:
        current_df.drop_duplicates(subset=primary_key, keep='last', inplace=True)
    else:
        current_df.drop_duplicates(keep='last', inplace=True)
    os.system(f'mv {path_to_current_file} {full_path}/backup/')
    os.system(f'cp "{path_to_new_table}" {full_path}/backup/tables/{os.path.basename(path_to_new_table)}')
    current_df.to_csv(f'{full_path}/{summary_name}_{report_time}{history_ext}', header=True, index=False)
    os.system(f'chmod 775 {full_path}/{summary_name}_{report_time}{history_ext}')


def update_files(
        arg_dict:dict, 
        parser, 
        current_tables:dict, 
        report_time:str, 
        full_path:str,
        proc_dict:dict,
        ignore_list:list=[]
        ):
    '''
    Wrapper function to perform history file update in a loop.
    arg_dict - refers to arg_dict global variable redefind in specific module;
    parser - argparse namespace object instantiated from the command line args;
    current_tables - dictionary mapping <report_prefix>_summary strings to tuple of form (path to current history file, dataframe of current history file))
    full_path - refers to full_path global variable redefind in specific module;
    report_time - refers to report_time global variable redefind in specific module (timestamp string expected);
    proc_dict - dictionary as it is defined in this module (proc_dict global variable), to supply case-specific settings to the process_report function (see README.md and modules for details)
    ignore_list - list of command-line arguments that do not refer to paths to file (e.g. batch number in update_aquamis_history.py);
    '''
    for report, arg in zip(arg_dict.keys(), vars(parser)):
        if report not in ignore_list:
            path_to_new_table    = getattr(parser,arg)
            path_to_current_file = current_tables[f'{report}_summary'][0]
            current_combined_df  = current_tables[f'{report}_summary'][1]
            summary_name         = f'{report}_summary'
            if arg not in proc_dict:
                pk = proc_dict['default']['primary_key']
                he = proc_dict['default']['history_ext']
                new_df = process_report(
                    path_to_new_table,
                    delimiter=proc_dict['default']['delimiter'],
                    col_rename_dict=proc_dict['default']['col_rename_dict'],
                    batch_regex=proc_dict['default']['batch_regex'],
                    id_column=proc_dict['default']['id_column'],
                    id_regex=proc_dict['default']['id_regex'],
                    batch_column="analysis_batch_id",
                    insert_batch=None,
                    custom_functions=proc_dict['default']['custom_functions']
                    )
            else:
                pk = proc_dict[arg]['primary_key']
                he = proc_dict[arg]['history_ext']
                new_df = process_report(
                    path_to_new_table,
                    delimiter=proc_dict[arg]['delimiter'],
                    col_rename_dict=proc_dict[arg]['col_rename_dict'],
                    batch_regex=proc_dict[arg]['batch_regex'],
                    id_column=proc_dict[arg]['id_column'],
                    id_regex=proc_dict[arg]['id_regex'],
                    batch_column=proc_dict[arg]['batch_column'],
                    insert_batch=proc_dict[arg]['insert_batch'],
                    custom_functions=proc_dict[arg]['custom_functions']
                    )
            update_history_file(
                path_to_new_table,
                path_to_current_file,
                current_combined_df,
                new_df,
                summary_name,
                report_time,
                full_path,
                primary_key=pk,
                history_ext=he
            )