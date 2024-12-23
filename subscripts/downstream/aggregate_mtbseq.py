import pandas as pd, time, os, sys, argparse, glob
from pathlib import Path
pd.set_option('mode.chained_assignment', None)

#################
#Global variables
#################

full_path = str(Path(os.path.dirname(os.path.realpath(__file__)))) #path to scripts
report_time = time.strftime("%m_%d_%Y") #timestamp
arg_dict = {
    "folder" : ["--fd" , "Folder containing MTBseq output folders. Each subfolder must include Classification, Called & Statistics subfolders."],
}

##########
#Functions
##########

def parse_arguments() -> argparse.ArgumentParser:
    '''Parsing predefined command-line arguments and return argparse class object.'''

    parser = argparse.ArgumentParser(description='A script to extract classification, statistics and variants from MTBseq output for single batch.') 
    
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


def process_report(path_to_df:str, delimiter:str='.tab', separator:str='\t') -> pd.DataFrame:
    """ 
    Returns cleaned pandas dataframe. 
    """

    if path_to_df.endswith(delimiter):
        try:
            df = pd.read_csv(path_to_df, sep=separator)
            return df
        except pd.errors.EmptyDataError:
            return pd.DataFrame()
    else:
        raise ValueError(f'{path_to_df} is expected to end with {delimiter}')


def clean_dataframe(df):
    '''
    Remove leading/trailing single quotes from string values
    '''
    for column in df.select_dtypes(include=['object']).columns:
        df[column] = df[column].str.strip("'")
    return df


def ET_mtbseq(file_map: dict, parser) -> None:
    '''
    Extract relevant information from MTBseq output files.
    '''
    # Initialize empty DataFrames for each report type
    stats_df = pd.DataFrame()
    class_df = pd.DataFrame()
    vars_df = pd.DataFrame()
    
    # Process statistics files
    for stats_file in file_map['stats']:
        batch_df = pd.read_csv(stats_file, sep='\t')
        batch_df = clean_dataframe(batch_df)
        stats_df = pd.concat([stats_df, batch_df], ignore_index=True)
    
    # Process classification files
    for class_file in file_map['class']:
        batch_df = pd.read_csv(class_file, sep='\t')
        batch_df = clean_dataframe(batch_df)
        class_df = pd.concat([class_df, batch_df], ignore_index=True)
    
    # Process variant files
    for vars_file in file_map['vars']:
        batch_df = pd.read_csv(vars_file, sep='\t')
        vars_df = pd.concat([vars_df, batch_df], ignore_index=True)
    
    # Create reports directory if it doesn't exist
    os.makedirs('reports', exist_ok=True)
    
    # Save combined reports
    stats_df.to_csv('reports/statistics_report.csv', index=False)
    class_df.to_csv('reports/classification_report.csv', index=False)
    vars_df.to_csv('reports/variants_report.csv', index=False)


def map_files(folder_path: str) -> dict:
    '''
    Returns dictionary mapping report type to full path as str: list.
    '''
    file_map = {
        'stats': os.path.join(folder_path, 'Statistics/Mapping_and_Variant_Statistics.tab'),
        'class': os.path.join(folder_path, 'Classification/Strain_Classification.tab'),
        'vars': os.path.join(folder_path, 'Called/*gatk_position_variants_cf4_cr4_fr75_ph4_outmode000.tab')
    }
    # Convert paths to lists of matching files
    return {k: [path for path in glob.glob(v)] for k, v in file_map.items()}

def get_mtbseq_folders(parent_folder: str) -> list:
    '''
    Returns list of valid MTBseq output folders that contain required subfolders.
    '''
    folders = []
    for folder in os.listdir(parent_folder):
        folder_path = os.path.join(parent_folder, folder)
        if os.path.isdir(folder_path):
            required = ['Statistics', 'Classification', 'Called']
            if all(os.path.isdir(os.path.join(folder_path, subfolder)) for subfolder in required):
                folders.append(folder_path)
    return folders




##############
#Runtime logic
##############

if __name__ == '__main__':
    parser = parse_arguments()
    
    # Get all valid MTBseq folders
    mtbseq_folders = get_mtbseq_folders(parser.fd)
    
    if not mtbseq_folders:
        print(f"No valid MTBseq output folders found in {parser.fd}")
        sys.exit(1)
    
    # Create combined file map with simpler keys
    combined_file_map = {
        'stats': [],
        'class': [],
        'vars': []
    }
    
    # Collect all files from each folder
    for folder in mtbseq_folders:
        folder_files = map_files(folder)
        for key in combined_file_map:
            combined_file_map[key].extend(folder_files[key])
    
    # Process all files at once
    ET_mtbseq(combined_file_map, parser)

