#!/opt/exp_soft/conda/anaconda3/bin/python

'''
This script is used to reorganize the contents of bact_output/*/ardetype/ folders.
Options for organizing the files:
a) By taxonomy
    - ardetype files for each sample will be placed according to its species (inferred from contings by kraken2)
    - aquamis files will stay unaltered in batch folders
b) By batch
    - all files are placed under corresponding batch subfolder
'''

###IMPORTS
import os, pandas as pd, sys, argparse, glob


###STATIC VARIABLES
ardetype_history_path = '/mnt/beegfs2/home/groups/nmrl/bact_analysis/analysis_history/ardetype/'
ardetype_file_tag = 'ardetype_history_file'
bact_output_path = '/mnt/beegfs2/home/groups/nmrl/bact_analysis/bact_output/'
taxonomy_parent_path = '/mnt/beegfs2/home/groups/nmrl/bact_analysis/bact_output/by_taxonomy/'


###HELPER FUNCTIONS
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    from https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()


def parse_arguments() -> argparse.ArgumentParser:
    '''Parsing cmd arguments and return argparse class object.'''

    #CMD ARGUMENTS & SCRIPT USAGE MESSAGES
    parser = argparse.ArgumentParser(description='A script structure ardetype output based on analysis_batch or taxonomy.\nWARNING: Make sure to update the summary report file!\nOnly data for samples present in summary file will be restructured (anything not present will stay in batch-based convention).') 
    parser.add_argument('-t', '--tax', help = 'Structure ardetype output by taxonomy (see by_taxonomy folder)', action='store_true')
    parser.add_argument('-b', '--batch', help = f'Structure ardetype output by batch (all folded outputs in corresponding batch folders)', action='store_true')


    #IF SCRIPT IS RUN WITHOUT ARGUMENTS - PRINT HELP MESSAGE & EXIT
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()


def find_history_file(folder_to_parse:str, tag_string:str) -> str:
    '''Parses the specified folder to find file path that contains tag string. Returns the string or None if not found.'''
    db_file_name = [file for file in os.listdir(folder_to_parse) if tag_string in file][0] #LOOKUP FOR THE DATABASE FILE IN THE FOLDER WHERE SCRIPT IS LOCATED
    if folder_to_parse[-1] == "/":
        db_file_path = folder_to_parse + db_file_name
    else:
        db_file_path = f'{folder_to_parse}/{db_file_name}'
    return db_file_path


def lower_unspace(string_list:list, repl_space:str="_", lowercase:bool=True) -> list:
    '''Given list of strings, replaces spaces and switches letters to lower case. Returns updated list.'''
    replaced_str = ";".join(string_list).replace(" ", repl_space)
    if lowercase: replaced_str = replaced_str.lower()
    return replaced_str.split(";")


def update_folders(name_list:list, subfolder_path:str, linux_perm:int=775) -> None:
    '''Given list of names, creates subfolders under <subfolder_path> directory, if folder with such name does not exist.'''
    for name in name_list: os.system(f'mkdir -m {linux_perm} -p {subfolder_path}{name}')


if __name__ == "__main__":
    curr_hist_file = find_history_file(ardetype_history_path, ardetype_file_tag)
    db_df = pd.read_csv(curr_hist_file)
    args = parse_arguments()

    if args.tax:
        folder_names = list(db_df['species'].unique())
        folder_names = lower_unspace(folder_names)
        update_folders(folder_names, taxonomy_parent_path)
        test_name = f"{taxonomy_parent_path}{db_df.loc[db_df['sample_id'] == db_df['sample_id'][0]]['species'].item().replace(' ', '_').lower()}/folded_{db_df['sample_id'][0]}_*_output/"
        if not glob.glob(test_name):
            print('Switching to taxonomy-based folder tree: ')
            for i,id in enumerate(db_df['sample_id']):
                source_name = f"{bact_output_path}{db_df.loc[db_df['sample_id'] == id]['analysis_batch_id'].item()}/*ardetype*/folded_{id}_*_output"
                dest_name = f"{taxonomy_parent_path}{db_df.loc[db_df['sample_id'] == id]['species'].item().replace(' ', '_').lower()}"
                os.system(f'mv {source_name} {dest_name} 2>/dev/null')
                printProgressBar(i+1, db_df['sample_id'].size)
        else:
            print('Already in taxonomy-based arrangement!')


    if args.batch:
        test_name = f"{taxonomy_parent_path}{db_df.loc[db_df['sample_id'] == db_df['sample_id'][0]]['species'].item().replace(' ', '_').lower()}/folded_{db_df['sample_id'][0]}_*_output/"
        if glob.glob(test_name):
            print('Switching to batch-based folder tree: ')
            for i,id in enumerate(db_df['sample_id']):
                dest_name = f"{bact_output_path}{db_df.loc[db_df['sample_id'] == id]['analysis_batch_id'].item()}/*ardetype*/"
                source_name = f"{taxonomy_parent_path}{db_df.loc[db_df['sample_id'] == id]['species'].item().replace(' ', '_').lower()}/folded_{id}_*_output"
                os.system(f'mv {source_name} {dest_name} 2>/dev/null')
                printProgressBar(i+1, db_df['sample_id'].size)
        else:
            print('Already in batch-based arrangement!')

