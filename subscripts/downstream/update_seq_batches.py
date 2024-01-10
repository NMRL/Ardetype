import sys, pandas as pd, time, os, pathlib, numpy as np
from datetime import datetime
sys.path.insert(0, f'/mnt/beegfs2/home/groups/nmrl/bact_analysis/Ardetype/subscripts/')
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

full_path   = str(pathlib.Path(os.path.dirname(os.path.realpath(__file__)))) #path to scripts



##########
#Functions
##########

def update_seq_batch() -> pd.DataFrame:
    '''
    Reads summary files from aquamis & ardetype subfolders.
    Combines unique identifiers (sample_id+analysis_batch_id) into separate tables.
    Pipeline columns, indicating the summary file of origin.
    '''
    #get paths to representative folder
    ardetype_path    = os.path.join(full_path,'ardetype')
    aquamis_path     = os.path.join(full_path,'aquamis')
    pointfinder_path = os.path.join(full_path,'resistance') #to include ids tagged with _reads
    mtbseq_path      = os.path.join(full_path,'mtbseq')     #to include mtbseq-only batches
    mobtyper_path    = os.path.join(full_path,'plasmids')     #to include mtbseq-only batches

    #get paths to current summaries
    ardetype_summary_path = next(pathlib.Path(ardetype_path).rglob('*ardetype_summary*.csv'), None)
    aquamis_summary_path  = next(pathlib.Path(aquamis_path).rglob('*aquamis_summary*.csv'), None) 
    pfr_summary_path      = next(pathlib.Path(pointfinder_path).rglob('*pointfinder_summary*.csv'), None)
    mtbseq_summary_path   = next(pathlib.Path(mtbseq_path).rglob('*classification_summary*.csv'), None)
    mobtyper_summary_path = next(pathlib.Path(mobtyper_path).rglob('*mobtyper_summary*.csv'), None)  

    #read summaries
    ar_df = pd.read_csv(ardetype_summary_path)
    aq_df = pd.read_csv(aquamis_summary_path)
    pf_df = pd.read_csv(pfr_summary_path)
    mb_df = pd.read_csv(mtbseq_summary_path)
    mt_df = pd.read_csv(mobtyper_summary_path)

    #preprocess summaries
    ar_b = ar_df[['sample_id', 'analysis_batch_id']]
    ar_b['tag'], ar_b['timestamp']  = [np.nan for _ in ar_b.index], [datetime.now().strftime('%Y-%m-%d-%H-%M') for _ in ar_b.index]
    aq_b = aq_df[['Sample_Name', 'analysis_batch_id']]
    aq_b = aq_b.rename(columns={'Sample_Name':'sample_id'})
    aq_b['tag'], aq_b['timestamp']  = [np.nan for _ in aq_b.index], [datetime.now().strftime('%Y-%m-%d-%H-%M') for _ in aq_b.index]
    pf_b = pf_df[['sample_id', 'analysis_batch_id']]
    pf_b = pf_b[pf_b.sample_id.str.contains('reads')]
    pf_b['tag'], pf_b['timestamp']  = [np.nan for _ in pf_b.index], [datetime.now().strftime('%Y-%m-%d-%H-%M') for _ in pf_b.index]
    mb_b = mb_df[['sample_id', 'analysis_batch_id']]
    mb_b['tag'], mb_b['timestamp']  = [np.nan for _ in mb_b.index], [datetime.now().strftime('%Y-%m-%d-%H-%M') for _ in mb_b.index]
    mt_df = mt_df[['sample_id', 'analysis_batch_id']]
    mt_df['tag'], mt_df['timestamp']  = [np.nan for _ in mt_df.index], [datetime.now().strftime('%Y-%m-%d-%H-%M') for _ in mt_df.index]

    #generate map
    seq_map = pd.read_csv(os.path.join(full_path,'seq_batch_map.csv'))
    seq_map = pd.concat([seq_map, ar_b, aq_b, pf_b, mt_df]).reset_index(drop=True).drop_duplicates(subset=['sample_id', 'analysis_batch_id'], keep='first').reset_index(drop=True)
    return seq_map


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
    summary_iterator = pathlib.Path('./').rglob('*_summary*.csv')
    to_exclude = ['backup', 'old', 'aquamis', 'ardetype', 'software']
    seq_batch_map['key'] = seq_batch_map['sample_id'].astype(str).str.cat(seq_batch_map['analysis_batch_id'].astype(str), sep='_')
    seq_batch_map_set = set(seq_batch_map['key'])

    summary_in_map_list = []
    map_in_summary_list = []

    for f in summary_iterator:
        if not any(s in str(f) for s in to_exclude):
            df = pd.read_csv(f, low_memory=False)
            df = df.dropna(subset=['sample_id', 'analysis_batch_id'])
            df['sample_id'] = df['sample_id'].astype(str)
            df['analysis_batch_id'] = df['analysis_batch_id'].astype(str)
            df['key'] = df['sample_id'].str.cat(df['analysis_batch_id'], sep='_')
            df_set = set(df['key'])

            missing_in_map = df_set - seq_batch_map_set
            missing_in_map_df = pd.DataFrame({'key': list(missing_in_map), 'present': str(f)})
            summary_in_map_list.append(missing_in_map_df)

            missing_in_summary = seq_batch_map_set - df_set
            missing_in_summary_df = pd.DataFrame({'key': list(missing_in_summary), 'missing': str(f)})
            map_in_summary_list.append(missing_in_summary_df)

    summary_in_map = pd.concat(summary_in_map_list).reset_index(drop=True)
    map_in_summary = pd.concat(map_in_summary_list).reset_index(drop=True)

    return summary_in_map, map_in_summary


def apply_tags(tag_file_path:str, sbm_path:str) -> pd.DataFrame:
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
    tag_file['timestamp'] = [datetime.now().strftime('%Y-%m-%d-%H-%M') for _ in tag_file.index]
    #combine tags and seq_batch_map
    merge = seq_batch_map.merge(tag_file, how='left', on=['sample_id', 'analysis_batch_id'])
    #Apply column merging logic described in the docstring
    merge['tag'] = merge['tag_y'].combine_first(merge['tag_x'])
    merge['timestamp'] = merge['timestamp_y'].combine_first(merge['timestamp_x'])
    #remove helper columns
    merge = merge.drop(['tag_x','tag_y', 'timestamp_x', 'timestamp_y'], axis=1)
    return merge

###############
# Runtime logic
###############

if __name__ == '__main__':
    args = hk.parse_arguments(arg_dict)
    if args.mode == 'update':
        new_map = update_seq_batch()
        new_map.to_csv(f'seq_batch_map.csv', header=True, index=False)
        sum_in_map, map_in_sum = check_batch_presence(new_map)
        sum_in_map.to_csv(f'missing_in_map.csv', header=True, index=False)
        map_in_sum.to_csv(f'missing_in_sum.csv', header=True, index=False)
    elif args.mode == 'tagging':
        #get full path to seq_batch_map.csv
        sbm_path = os.path.join(full_path,'seq_batch_map.csv')
        #apply tags from file
        tagged_sbm = apply_tags(tag_file_path=args.tag_file, sbm_path=sbm_path)
        #save tagged seq_batch_map to the original path
        tagged_sbm.to_csv(sbm_path, header=True, index=False)