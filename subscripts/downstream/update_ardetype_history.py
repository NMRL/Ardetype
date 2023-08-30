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
        help = 'Full path to ardetype report (e.g. ardetype_report.csv)', 
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


def process_report(ardetype_report_path:str) -> pd.DataFrame:
    '''
    Given path to ardetype report, restructures it by concatenating method, 
    type and reference columns for each sample into single column, sparated by doubly-spaced semicolon ( ; ).
    If all fields are empty for any given sample, empty columns are produced for that sample. 
    Removed illumina sample number from all sample names.
    Returns reformatted pandas dataframe.
    '''
    df = pd.read_csv(ardetype_report_path).applymap(str)
    df['sample_id'] = df['sample_id'].str.replace(r"_S[0-9]*","")
    tools = [col.replace("method|", "") for col in df.columns if "method" in col]
    if tools: #IF SPECIES-SPECIFIC TYPING WAS PERFORMED
        df["methods"] = df[[col for col in df.columns if "method|" in col]].agg(' ; '.join, axis=1) #CONCATENATE METHODS
        df["types"] = df[[col for col in df.columns if "type|" in col]].agg(' ; '.join, axis=1) #CONCATENATE TYPES
        df["types"] = [(' ; ').join([tools[i]+": "+value for i,value in enumerate(row.split(' ; '))]) for row in df['types']] #ADD TOOL-BASED ANNOTATION TO EACH TYPE
        df["references"] = df[[col for col in df.columns if "reference|" in col]].agg(' ; '.join, axis=1) #CONCATENATE REFERENCES
        df.drop(list(df.filter(regex = '(method\||type\||reference\|)')), axis = 1, inplace = True) #REMOVE UNUSED COLUMNS
        df = df.replace('((; |)nan( |)|; [A-z0-9]*: nan |[A-z0-9]*: nan ; |[\-A-z0-9]*: nan)','', regex=True) #REMOVE REDUNDANT INFORMATION
        df = df.replace('(; $|^; )','', regex=True) #REMOVE REDUNDANT SEPARATORS
        df.drop('taxid', inplace=True, axis=1) #DROP UNUSED COLUMN
    else:
        #BLANK ENTRIES FOR SAMPLES WHERE NO SPECIFIC DATA IS AVAILABLE
        df['methods'] = ["" for _ in df.index]
        df['types'] = ["" for _ in df.index] 
        df['references'] = ["" for _ in df.index]
        df.fillna("", inplace=True)
    return df


def find_current_table():
    '''
    Checks files in the folder where the script is placed. 
    Returns tuple containing pd.Dataframe object generated from that file and path to the file, if file name contains "ardetype_summary" substring. 
    Returns None if file is not found.
    '''
    ardetype_summary, adt_df      = '', None
    for file in os.listdir(full_path):
        path_to_file = f'{full_path}/{file}'
        if os.path.isfile(path_to_file) and 'ardetype_summary' in file: 
            ardetype_summary, adt_df      = path_to_file, pd.read_csv(path_to_file, low_memory=False)
    if adt_df is not None:
        return ardetype_summary, adt_df
    else:
        return None


def update_ardetype_file(path_to_new_table:str, path_to_current_file:str, current_combined_df:pd.DataFrame, new_df:pd.DataFrame, type='ardetype'):
    '''
    Given dataframe containing formatted batch of ardetype, 
    adds it to the combined dataframe, 
    moves previously generated combined table file to backup/, 
    moves new resistance table to backyp/weekly/,
    saves updated combined dataframe as "ardetype_summary" file with current timestamp and sets its permissions to 775.
    '''
    if type == 'ardetype':
        current_df = pd.concat([current_combined_df,new_df], sort=False)
        current_df = current_df.astype(str)
        current_df.drop_duplicates(subset=['sample_id', 'analysis_batch_id'], keep='last', inplace=True)
        os.system(f'mv {path_to_current_file} {full_path}/backup/')
        os.system(f'cp "{path_to_new_table}" {full_path}/backup/tables/{os.path.basename(path_to_new_table)}')
        current_df.to_csv(f'{full_path}/ardetype_summary_{report_time}.csv', header=True, index=False)
        os.system(f'chmod 775 {full_path}/ardetype_summary_{report_time}.csv')


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
    adt_df = process_report(os.path.abspath(parser.pipeline)) 
    #update hamronization-aggregated summary table
    update_ardetype_file( 
        path_to_new_table=parser.pipeline ,
        current_combined_df=current_combined_tuple[1], 
        path_to_current_file=current_combined_tuple[0],
        new_df=adt_df
    )
