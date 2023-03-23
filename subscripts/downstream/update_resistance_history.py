import pandas as pd, time, os, sys
from pathlib import Path
pd.set_option('mode.chained_assignment', None) 


def create_backup():
    '''Create backup folder if it does not exist'''
    if not os.path.isdir(f'{full_path}/backup/tables'):
        os.system(f'mkdir -p {full_path}/backup/tables')


def process_report(path_to_df:str) -> pd.DataFrame:
    """ Reads resistance report into pandas dataframe. Removes illumina sample number from sample ids. Returns cleaned dataframe. """
    if path_to_df.endswith(".tsv"):
        df = pd.read_csv(path_to_df, sep='\t')
        df.input_file_name = df.input_file_name.str.replace('.amr.alignment', '')
        df.input_file_name = df.input_file_name.str.replace(r"(_S[0-9]{3}|_S[0-9]{2}|_S[0-9]{1})", '')
        return df
    elif path_to_df.endswith('.html'):
        pass

def find_current_table(base_path:str):
    '''
    Checks files in the folder where the script is placed. 
    Returns tuple containing pd.Dataframe object generated from that file and path to the file, if file name contains "resistance_summary" substring. 
    Returns None if file is not found.
    '''
    resistance_summary, rs_df      = '', None
    resfinder_pheno_summar, rfp_df = '', None
    for file in os.listdir(base_path):
        path_to_file = f'{base_path}/{file}'
        if os.path.isfile(path_to_file) and 'resistance_summary' in file: 
            resistance_summary, rs_df      = path_to_file, pd.read_csv(file, low_memory=False, sep='\t')
        elif os.path.isfile(path_to_file) and 'resfinder_pheno_summary' in file:
            resfinder_pheno_summar, rfp_df = path_to_file, pd.read_csv(file, low_memory=False)
    if rs_df is not None and rfp_df is not None:
        return resistance_summary, rs_df, resfinder_pheno_summar, rfp_df
    else:
        return None


def update_combined_resistance_file(base_path:str, report_time:str, path_to_new_table:str, path_to_current_file:str, current_combined_df:pd.DataFrame, new_df:pd.DataFrame, type='resistance'):
    '''
    Given dataframe containing formatted batch of resistance, adds it to the combined dataframe, 
    moves previously generated combined table file to backup/, 
    moves new resistance table to backyp/weekly/,
    saves updated combined dataframe as "resistance_summary" file with current timestamp and sets its permissions to 775.
    '''
    if type == 'resistance':
        current_df = pd.concat([current_combined_df,new_df], sort=False)
        current_df.drop_duplicates(keep='first', inplace=True)
        os.system(f'mv {path_to_current_file} {base_path}/backup/')
        os.system(f'cp "{path_to_new_table}" {base_path}/backup/tables/{os.path.basename(path_to_new_table).replace("(sequencing_date)_(batch)",time.strftime("%Y_%m_%d:%H:%M:%S"))}')
        current_df.to_csv(f'{base_path}/resistance_summary_{report_time}.tsv', sep='\t', header=True, index=False)
        os.system(f'chmod 775 {base_path}/resistance_summary_{report_time}.tsv')
    elif type == 'resfinder':
        current_df = pd.concat([current_combined_df,new_df], sort=False)
        current_df.drop_duplicates(keep='first', inplace=True)
        os.system(f'mv {path_to_current_file} {base_path}/backup/')
        os.system(f'cp "{path_to_new_table}" {base_path}/backup/tables/{os.path.basename(path_to_new_table).replace(".csv",time.strftime("%Y_%m_%d:%H:%M:%S.csv"))}')
        current_df.to_csv(f'{base_path}/resfinder_summary_{report_time}.csv', header=True, index=False)
        os.system(f'chmod 775 {base_path}/resfinder_summary_{report_time}.csv')


if __name__ == '__main__':
    full_path = str(Path(os.path.dirname(os.path.realpath(__file__)))) #path to scripts
    report_time = time.strftime("%m_%d_%Y") #timestamp
    args = sys.argv #positional cmd arguments
    create_backup() #if not present
    current_combined_tuple = find_current_table(full_path) #find and read current summary table
    if len(args) > 1 is not None:
        res_df = process_report(args[1]) #find, read and preprocess new table
        rfp_df = pd.read_csv(args[2])
        update_combined_resistance_file( #update resistance
            base_path=full_path,
            report_time=report_time,
            path_to_new_table=args[1], 
            current_combined_df=current_combined_tuple[1], 
            path_to_current_file=current_combined_tuple[0],
            new_df=res_df
        )
        update_combined_resistance_file(
            base_path=full_path,
            report_time=report_time,
            path_to_new_table=args[2], 
            path_to_current_file=current_combined_tuple[2], 
            current_combined_df=current_combined_tuple[3], 
            new_df=rfp_df,
            type="resfinder"
        )
    else:
        sys.exit(f'''
        Please ensure that both resfinder_pheno_summary file and resistance_summary file are present.
        ''')


