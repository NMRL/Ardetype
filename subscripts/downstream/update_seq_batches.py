import pandas as pd, time, os, pathlib
pd.set_option('mode.chained_assignment', None)

#################
#Global variables
#################

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

    #get paths to current summaries
    ardetype_summary_path = [f for f in pathlib.Path(ardetype_path).rglob('*ardetype_summary*.csv')][0]
    aquamis_summary_path  = [f for f in pathlib.Path(aquamis_path).rglob('*aquamis_summary*.csv')][0]
    pfr_summary_path      = [f for f in pathlib.Path(pointfinder_path).rglob('*pointfinder_summary*.csv')][0]
    mtbseq_summary_path   = [f for f in pathlib.Path(mtbseq_path).rglob('*classification_summary*.csv')][0]

    #read summaries
    ar_df = pd.read_csv(ardetype_summary_path)
    aq_df = pd.read_csv(aquamis_summary_path)
    pf_df = pd.read_csv(pfr_summary_path)
    mb_df = pd.read_csv(mtbseq_summary_path)

    #preprocess summaries
    ar_b = ar_df[['sample_id', 'analysis_batch_id']]
    aq_b = aq_df[['Sample_Name', 'analysis_batch_id']] 
    aq_b = aq_b.rename(columns={'Sample_Name':'sample_id'})
    pf_b = pf_df[['sample_id', 'analysis_batch_id']]
    pf_b = pf_b[pf_b.sample_id.str.contains('reads')]
    mb_b = mb_df[['sample_id', 'analysis_batch_id']]

    #generate map
    seq_map = pd.concat([ar_b, aq_b, pf_b]).reset_index(drop=True).drop_duplicates(keep='last').reset_index(drop=True)
    return seq_map


def check_batch_presence(seq_batch_map:pd.DataFrame) -> tuple:
    '''
    Lists all relevant summary files.
    Parses as dataframes.
    For each batch id (sample_id+analysis_batch_id) in seq_batch_map, checks its presence in corresponding dataframe.
    For each batch id in every dataframe, checks its presence in seq_batch_map.
    Returns two dataframes.
    First contains list of unique batch ids from seq_batch_map and corresponding status in every summary dataframe.
    Second contains list of unique batch ids from every dataframe and corresponding status in seq_batch_map.
    '''
    summary_iterator     = pathlib.Path('./').rglob('*_summary*.csv')
    to_exclude           = ['backup', 'old', 'aquamis', 'ardetype', 'software']
    seq_batch_map['key'] = seq_batch_map.sample_id.str.cat(seq_batch_map.analysis_batch_id, sep='_')

    summary_in_map = pd.DataFrame()
    map_in_summary = pd.DataFrame()
    for f in summary_iterator:
        check = any([s in str(f) for s in to_exclude])
        if not check:
            #preprocess the file
            path                      = str(f)
            df                        = pd.read_csv(path, low_memory=False)
            df                        = df.dropna(subset=['sample_id', 'analysis_batch_id'])
            df                        = df.astype(str)
            df_b                      = df[['sample_id', 'analysis_batch_id']]
            df_b                      = df_b.drop_duplicates(keep='last')
            df_b['key']               = df_b.sample_id.str.cat(df_b.analysis_batch_id, sep='_')
            #find missing in summaries
            missing_in_map            = df_b[~df_b.key.isin(seq_batch_map.key)].reset_index(drop=True)
            missing_in_map['present'] = [path for _ in missing_in_map.index]
            #find missing in seq_map
            missing_in_sum            = seq_batch_map[~seq_batch_map.key.isin(df_b.key)].drop(['key'],axis=1).reset_index(drop=True)
            missing_in_sum['missing'] = [path for _ in missing_in_sum.index]
            #aggragate
            summary_in_map            = pd.concat([summary_in_map, missing_in_map])
            map_in_summary            = pd.concat([map_in_summary, missing_in_sum])

    return summary_in_map, map_in_summary


# ##############
# #Runtime logic
# ##############

if __name__ == '__main__':

    new_map = update_seq_batch()
    new_map.to_csv(f'seq_batch_map.csv', header=True, index=False)
    sum_in_map, map_in_sum = check_batch_presence(new_map)
    sum_in_map.to_csv(f'missing_in_map.csv', header=True, index=False)
    map_in_sum.to_csv(f'missing_in_sum.csv', header=True, index=False)

