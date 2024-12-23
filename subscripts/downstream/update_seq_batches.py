import sys, pandas as pd, os, pathlib
import concurrent.futures
from datetime import datetime
sys.path.insert(0, f'/home/group/pipelines/Ardetype/subscripts/')
from src.utilities import Housekeeper as hk
pd.set_option('mode.chained_assignment', None)


#################
#Global variables
#################

#defining cmd arguments
arg_dict = {
    "description": "",
    "required_arguments": [
        [
            "-m",
            "--mode",
            "Selecting mode to run the script:\nupdate - scan folders and update list of identifiers \ntagging - apply tags from external csv files."
        ]
    ],
    "optional": {
        "arguments": [
            [
                "-tf",
                "--tag_file",
                "Path to a csv file that maps tags to sample identifier (sample_id, analysis_batch_id) - requires <tagging> mode.",
                None
            ]
        ],
        "flags": [
            # [
            #     "-s",
            #     "--install_snakemake",
            #     "Flag to install mamba and snakemake for the current HPC user, if it is not already installed."
            # ]
        ],
        "nargs": [
            # [
            #     "-mf",
            #     "--merge_from",
            #     "Argument to supply paths to folders containing ardetype output that should be combined."
            # ]
        ]
    }
}

#global static variables
full_path   = str(pathlib.Path(os.path.dirname(os.path.realpath(__file__)))) #path to scripts
date = datetime.now().isoformat(' ', 'milliseconds')  #1st argument defines date and time separator, 2nd - precision - see https://www.geeksforgeeks.org/isoformat-method-of-datetime-class-in-python/

##########
#Functions
##########

def update_seq_batch(date:str, null_plh='None') -> pd.DataFrame:
    '''
    Reads summary files from aquamis & ardetype subfolders.
    Combines unique identifiers (sample_id+analysis_batch_id) into separate tables.
    Pipeline columns, indicating the summary file of origin.
    '''
    #name folders under analysis_history/ that contain reports for all samples in all batches
    #folder name is mapped to the name of the report file
    rep_map = {
        'ardetype':'ardetype',
        'resistance':'pointfinder', #to include ids tagged with _reads
        'mtbseq':'classification', #to include mtbseq-only batches
        'plasmids':'mobtyper' #to include mtbseq-only batches
    }
    #convert to full paths
    rep_fold_paths = [os.path.join(full_path, folder_name) for folder_name in rep_map]

    #get paths to current summaries
    cur_summary_paths = [next(pathlib.Path(rep_fold_path).rglob(f'*{rep_map[folder_name]}_summary*.csv')) for rep_fold_path, folder_name in zip(rep_fold_paths, rep_map.keys())]

    #read summaries
    summary_list = [pd.read_csv(path) for path in cur_summary_paths]
    
    #preprocess summaries
    def extract_batch_timestamp(
            df:pd.DataFrame, 
            batch_timestamp_pattern:str=r'(_[0-9]{8}_[0-9]{6})', 
            sbt__in_format:str='%Y%m%d_%H%M%S', 
            sbt_out_format:str='%Y-%m-%d_%H-%M-%S',
            null_placeholder:str=null_plh,
            ) -> pd.DataFrame:
        '''
        Extracting timestamp substring from seq_batch_id
        '''
        df['seq_batch_timestamp'] = df['analysis_batch_id'].str.extract(batch_timestamp_pattern)[0].str.lstrip('_')
        df['seq_batch_timestamp'] = pd.to_datetime(df['seq_batch_timestamp'], format=sbt__in_format).dt.strftime(sbt_out_format).astype(str)
        df['seq_batch_timestamp'] = df['seq_batch_timestamp'].replace('NaT', null_placeholder)
        return df

    #define column names to be normalized to single format
    special_id_cases = ['Sample_Name', 'input_file_name']
    special_batch_cases = []
    for i,df in enumerate(summary_list):
        #normalize column names
        for case in special_id_cases: df.rename(columns={case:'sample_id'}, inplace=True)
        for case in special_batch_cases: df.rename(columns={case:'analysis_batch_id'}, inplace=True)

        #keep GUIDs only
        df = df[['sample_id', 'analysis_batch_id']]

        #add tag and tag_timestamp for new records
        df['tag'], df['tag_timestamp']  = [null_plh for _ in df.index], [date for _ in df.index]

        #separate seq_batch_timestamp into separate column in defined format
        df = extract_batch_timestamp(df)
        summary_list[i] = df

    #get existing map
    seq_map = pd.read_csv(os.path.join(full_path,'seq_batch_map.csv'))
    #adding existing records to avoid overwriting timestamps
    summary_list.insert(0,seq_map)
    #combine existing map with any new records keeping existing records when duplicates are dropped
    seq_map = pd.concat(summary_list, sort=False).reset_index(drop=True).drop_duplicates(subset=['sample_id', 'analysis_batch_id'], keep='first').reset_index(drop=True)
    return seq_map

def process_file_wrapper(args):
    '''To make the function suitable for concurrency'''
    def process_file(f, seq_batch_map_set, sample_id_columns, to_exclude):
        '''GUID extractor to check presence in seq_batch_map'''
        if not any(s in str(f) for s in to_exclude):
            df = pd.read_csv(f, low_memory=False)
            #normalize id column to sample_id
            for c in sample_id_columns:
                if c in df.columns:
                    df.rename(columns={c:'sample_id'}, inplace=True)
            #remove records with missing GUIDS
            df = df.dropna(subset=['sample_id', 'analysis_batch_id'])
            #cast GUIDs to strings to apply cat operation
            df[['sample_id', 'analysis_batch_id']] = df[['sample_id', 'analysis_batch_id']].astype(str)
            #merge two-column GUID to single column
            df['key'] = df['sample_id'].str.cat(df['analysis_batch_id'], sep='_')
            #casting to set 
            df_set = set(df['key'])
            #detect ids present in summary but missing in map using set difference operation
            missing_in_map = df_set - seq_batch_map_set
            #keeping the results (f stands for summary file name)
            missing_in_map_df = pd.DataFrame({'key': list(missing_in_map), 'present': str(f)})
            
            #same for present in seq_batch_map but missing in specific report
            missing_in_summary = seq_batch_map_set - df_set
            missing_in_summary_df = pd.DataFrame({'key': list(missing_in_summary), 'missing': str(f)})
            
            return missing_in_map_df, missing_in_summary_df
        return pd.DataFrame(), pd.DataFrame()
    return process_file(*args)


def check_batch_presence(seq_batch_map: pd.DataFrame) -> tuple:
    '''
    Lists all relevant summary files.
    Parses as dataframes.
    For each batch id (sample_id+analysis_batch_id) in seq_batch_map, checks its presence in corresponding dataframe.
    For each batch id in every dataframe, checks its presence in seq_batch_map.
    Returns two dataframes.
    First contains list of unique batch ids from seq_batch_map and corresponding status in every summary dataframe.
    Second contains list of unique batch ids from every dataframe and corresponding status in seq_batch_map.
    '''
    summary_iterator = list(pathlib.Path('./').rglob('*_summary*.csv'))
    to_exclude = ['backup', 'old', 'aquamis', 'ardetype', 'software', 'k2reads']
    seq_batch_map['key'] = seq_batch_map['sample_id'].astype(str).str.cat(seq_batch_map['analysis_batch_id'].astype(str), sep='_')
    seq_batch_map_set = set(seq_batch_map['key'])
    sample_id_columns = ['input_file_name']

    summary_in_map_list = []
    map_in_summary_list = []

    with concurrent.futures.ProcessPoolExecutor() as executor:
        # Prepare the arguments for the helper function
        args = [(f, seq_batch_map_set, sample_id_columns, to_exclude) for f in summary_iterator]

        # Use the helper function with executor.map
        results = executor.map(process_file_wrapper, args)
        
        #aggregate per-summary results
        for missing_in_map_df, missing_in_summary_df in results:
            summary_in_map_list.append(missing_in_map_df)
            map_in_summary_list.append(missing_in_summary_df)

    #combine into single dataframe
    summary_in_map = pd.concat(summary_in_map_list).reset_index(drop=True)
    map_in_summary = pd.concat(map_in_summary_list).reset_index(drop=True)

    return summary_in_map, map_in_summary


def apply_tags(tag_file_path:str, sbm_path:str, date:str) -> pd.DataFrame:
    '''
    input:
        tag_file_path - full path to a csv file with the following structure
            sample_id,analysis_batch_id,tag
            ...
            <sample_id + analysis_batch_id> - unique record identifier
            <tag> - any tag to be recognized by the downstream applications
        sbm_path - full path to current seq_batch_map.csv file containing all identifiers currently present in all history files
        Empty tags format - empty cell in a csv file
    output:
        tagged seq_batch_map as pandas dataframe

    tagging logic
    Tag supplied in tag_file will have prevalence over existing tag (overwrites it) !if it is not an empty tag!.
    Examples of tagging:
    1)
    original tag: empty (no value)
    incoming tag: string A
    resulting tag: string A
    2) 
    original tag: string A
    incoming tag: string B
    resulting tag: string B
    3)
    original tag: string A
    incoming tag: empty (no value)
    resulting tag: string A
    4)
    original tag: empty (no value)
    incoming tag: empty (no value)
    resulting tag: empty (no value)
    '''

    #read with identifiers as strings for smooth merging
    tag_file = pd.read_csv(tag_file_path, dtype={'sample_id':str,'analysis_batch_id':str})
    seq_batch_map = pd.read_csv(sbm_path, dtype={'sample_id':str,'analysis_batch_id':str})
    #generate timestamps at the moment when the tags are applied
    tag_file['tag_timestamp'] = [date for _ in tag_file.index]
    #combine tags and seq_batch_map
    merge = seq_batch_map.merge(tag_file, how='left', on=['sample_id', 'analysis_batch_id'])
    #Apply column merging logic described in the docstring
    merge['tag'] = merge['tag_y'].combine_first(merge['tag_x'])
    merge['tag_timestamp'] = merge['tag_timestamp_y'].combine_first(merge['tag_timestamp_x'])
    #remove helper columns
    merge = merge.drop(['tag_x','tag_y', 'tag_timestamp_x', 'tag_timestamp_y'], axis=1)
    return merge

###############
# Runtime logic
###############

def main():
    args = hk.parse_arguments(arg_dict)
    if args.mode == 'update':
        new_map = update_seq_batch(date=date)
        new_map.to_csv(f'seq_batch_map.csv', header=True, index=False)
        sum_in_map, map_in_sum = check_batch_presence(new_map)
        sum_in_map.to_csv(f'missing_in_map.csv', header=True, index=False)
        map_in_sum.to_csv(f'missing_in_sum.csv', header=True, index=False)
    elif args.mode == 'tagging':
        #get full path to seq_batch_map.csv
        sbm_path = os.path.join(full_path,'seq_batch_map.csv')
        #apply tags from file
        tagged_sbm = apply_tags(tag_file_path=args.tag_file, sbm_path=sbm_path, date=date)
        #save tagged seq_batch_map to the original path
        tagged_sbm.to_csv(sbm_path, header=True, index=False)


if __name__ == '__main__':
    main()