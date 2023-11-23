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
        '-plf', 
        '--plasmidfinder', 
        metavar='\b', 
        help = 'Full path to the plasmidfinder summary report (e.g. plasmidfinder_summary.csv)', 
        default=None, 
        required=True
        )
    parser.add_argument(
        '-mbt', 
        '--mob_typer', 
        metavar='\b', 
        help = f'Full path to the aggregated mob_typer summary report (e.g. mobtyper_summary.csv)', 
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


def process_report(path_to_df:str, remove_sn:bool=False) -> pd.DataFrame:
    """ 
    Returns cleaned dataframe. 
    """

    if path_to_df.endswith(".csv"):
        df                 = pd.read_csv(path_to_df)
        df.sample_id = df.sample_id.astype(str)
        if remove_sn: df.sample_id = df.sample_id.str.replace(r"_S[0-9]*", '')
        df.analysis_batch_id = df.analysis_batch_id.str.replace(r'_[0-9]{8}_[0-9]{6}', '')
        return df
    else:
        raise ValueError(f'{path_to_df} is expected to end with .csv')


def find_current_table():
    '''
    Checks files in the folder where the script is placed. 
    Returns tuple containing pd.Dataframe object generated from that file and path to the file, 
    if file name contains "{plasmid_tool}_summary" substring. 
    Returns None if file is not found.
    '''
    plasmidfinder_summary, pf_df      = '', None
    mobtyper_summary, mt_df = '', None
    for file in os.listdir(full_path):
        path_to_file = f'{full_path}/{file}'
        if os.path.isfile(path_to_file) and 'plasmidfinder_summary' in file: 
            plasmidfinder_summary = path_to_file
            try:
                pf_df = pd.read_csv(path_to_file, low_memory=False)
            except pd.errors.EmptyDataError:
                pf_df = pd.DataFrame()
        elif os.path.isfile(path_to_file) and 'mobtyper_summary' in file:
            mobtyper_summary = path_to_file
            try:
                mt_df = pd.read_csv(path_to_file, low_memory=False)
            except pd.errors.EmptyDataError:
                mt_df = pd.DataFrame()
    if pf_df is not None and mt_df is not None:
        return plasmidfinder_summary, pf_df, mobtyper_summary, mt_df
    else:
        return None


def update_combined_plasmid_file(
        path_to_new_table:str, 
        path_to_current_file:str, 
        current_combined_df:pd.DataFrame, 
        new_df:pd.DataFrame, 
        type:str) -> None:
    '''
    Given dataframe containing formatted batch of plasmid, 
    adds it to the combined dataframe, 
    moves previously generated combined table file to backup/, 
    moves new plasmid table to backyp/weekly/,
    saves updated combined dataframe as "plasmid_summary" file with current timestamp and sets its permissions to 775.
    '''
    if type == 'plasmidfinder':
        current_df = pd.concat([current_combined_df,new_df], sort=False)
        current_df.drop_duplicates(keep='last', inplace=True)
        os.system(f'mv {path_to_current_file} {full_path}/backup/')
        os.system(f'cp "{path_to_new_table}" {full_path}/backup/tables/{os.path.basename(path_to_new_table).replace(".csv",time.strftime("%Y_%m_%d:%H:%M:%S"))}')
        current_df.to_csv(f'{full_path}/plasmidfinder_summary_{report_time}.csv', header=True, index=False)
        os.system(f'chmod 775 {full_path}/plasmidfinder_summary_{report_time}.csv')
    elif type == 'mobtyper':
        #separate contam.records
        current_df = pd.concat([current_combined_df,new_df], sort=False)
        current_df.drop_duplicates(keep='last', inplace=True)
        os.system(f'mv {path_to_current_file} {full_path}/backup/')
        os.system(f'cp "{path_to_new_table}" {full_path}/backup/tables/{os.path.basename(path_to_new_table).replace(".csv",time.strftime("_%Y_%m_%d:%H:%M:%S.csv"))}')
        current_df.to_csv(f'{full_path}/mobtyper_summary_{report_time}.csv', header=True, index=False)
        os.system(f'chmod 775 {full_path}/mobtyper_summary_{report_time}.csv')


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
    res_df = process_report(os.path.abspath(parser.plasmidfinder), remove_sn=True) 
    rfp_df = process_report(os.path.abspath(parser.mob_typer))
    
    #update plasmidfinder summary table
    update_combined_plasmid_file( 
        path_to_new_table=parser.plasmidfinder ,
        current_combined_df=current_combined_tuple[1], 
        path_to_current_file=current_combined_tuple[0],
        new_df=res_df,
        type='plasmidfinder'
    )

    #update mob_typer summary table
    update_combined_plasmid_file(
        path_to_new_table=parser.mob_typer, 
        path_to_current_file=current_combined_tuple[2], 
        current_combined_df=current_combined_tuple[3], 
        new_df=rfp_df,
        type="mobtyper"
    )
