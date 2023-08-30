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
        '-hamr', 
        '--hamronization', 
        metavar='\b', 
        help = 'Full path to the hamronization summary report (e.g. summarized_resistance_profile_(sequencing_date)_(batch).tsv)', 
        default=None, 
        required=True
        )
    parser.add_argument(
        '-rf', 
        '--resfinder', 
        metavar='\b', 
        help = f'Full path to the aggregated resfinder phenotype report (e.g. resfinder_pheno_table_gathered.csv)', 
        default=None, 
        required=True
        )
    parser.add_argument(
        '-pf', 
        '--pointfinder', 
        metavar='\b', 
        help = f'Full path to the aggregated pointfinder report (e.g. pointfinder_report.csv)', 
        default=None, 
        required=True
        )

    parser.add_argument(
        '-afm', 
        '--amrfinder_mut', 
        metavar='\b', 
        help = f'Full path to the aggregated amrfinder mutation report (e.g. amrfp_mutation_report.csv)', 
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


def process_report(path_to_df:str) -> pd.DataFrame:
    """ 
    Reads hamronization resistance report (in tsv format) into pandas dataframe. 
    Removes illumina sample number from sample ids. 
    Returns cleaned dataframe. 
    """

    if path_to_df.endswith(".tsv"):
        df                 = pd.read_csv(path_to_df, sep='\t')
        df.input_file_name = df.input_file_name.astype(str)
        df.input_file_name = df.input_file_name.str.replace('.amr.alignment', '')
        df.input_file_name = df.input_file_name.str.replace(r"_S[0-9]*", '')
        return df
    else:
        raise ValueError(f'{path_to_df} is expected to end with .tsv')


def find_current_table():
    '''
    Checks files in the folder where the script is placed. 
    Returns tuple containing pd.Dataframe object generated from that file and path to the file, if file name contains "resistance_summary" substring. 
    Returns None if file is not found.
    '''
    resistance_summary, rs_df      = '', None
    resfinder_pheno_summar, rfp_df = '', None
    for file in os.listdir(full_path):
        path_to_file = f'{full_path}/{file}'
        if os.path.isfile(path_to_file) and 'resistance_summary' in file: 
            resistance_summary, rs_df      = path_to_file, pd.read_csv(path_to_file, low_memory=False, sep='\t')
        elif os.path.isfile(path_to_file) and 'resfinder_summary' in file:
            resfinder_pheno_summar, rfp_df = path_to_file, pd.read_csv(path_to_file, low_memory=False)
        elif os.path.isfile(path_to_file) and 'pointfinder_summary' in file:
            pointfinder_summary, pf_df     = path_to_file, pd.read_csv(path_to_file, low_memory=False)
        elif os.path.isfile(path_to_file) and 'amrfinder_mut_summary' in file:
            amrfm_summary, afm_df     = path_to_file, pd.read_csv(path_to_file, low_memory=False)

    
    if rs_df is not None and rfp_df is not None and pf_df is not None and afm_df is not None:
        return resistance_summary, rs_df, resfinder_pheno_summar, rfp_df, pointfinder_summary, pf_df, amrfm_summary, afm_df
    else:
        return None


def update_combined_resistance_file(path_to_new_table:str, path_to_current_file:str, current_combined_df:pd.DataFrame, new_df:pd.DataFrame, type='resistance'):
    '''
    Given dataframe containing formatted batch of resistance, 
    adds it to the combined dataframe, 
    moves previously generated combined table file to backup/, 
    moves new resistance table to backyp/weekly/,
    saves updated combined dataframe as "resistance_summary" file with current timestamp and sets its permissions to 775.
    '''
    if type == 'resistance':
        current_df = pd.concat([current_combined_df,new_df], sort=False)
        current_df.drop_duplicates(keep='last', inplace=True)
        os.system(f'mv {path_to_current_file} {full_path}/backup/')
        os.system(f'cp "{path_to_new_table}" {full_path}/backup/tables/{os.path.basename(path_to_new_table).replace("(sequencing_date)_(batch)",time.strftime("%Y_%m_%d:%H:%M:%S"))}')
        current_df.to_csv(f'{full_path}/resistance_summary_{report_time}.tsv', sep='\t', header=True, index=False)
        os.system(f'chmod 775 {full_path}/resistance_summary_{report_time}.tsv')
    elif type == 'resfinder':
        #separate contam.records
        current_df = pd.concat([current_combined_df,new_df], sort=False)
        current_df.drop_duplicates(keep='last', inplace=True)
        current_df.dropna(subset=['Genetic background'], inplace=True)
        os.system(f'mv {path_to_current_file} {full_path}/backup/')
        os.system(f'cp "{path_to_new_table}" {full_path}/backup/tables/{os.path.basename(path_to_new_table).replace(".csv",time.strftime("_%Y_%m_%d:%H:%M:%S.csv"))}')
        current_df.to_csv(f'{full_path}/resfinder_summary_{report_time}.csv', header=True, index=False)
        os.system(f'chmod 775 {full_path}/resfinder_summary_{report_time}.csv')
    elif type == 'pointfinder':
        current_df = pd.concat([current_combined_df,new_df], sort=False)
        current_df.drop_duplicates(keep='last', inplace=True)
        os.system(f'mv {path_to_current_file} {full_path}/backup/')
        os.system(f'cp "{path_to_new_table}" {full_path}/backup/tables/{os.path.basename(path_to_new_table).replace(".csv",time.strftime("_%Y_%m_%d:%H:%M:%S.csv"))}')
        current_df.to_csv(f'{full_path}/pointfinder_summary_{report_time}.csv', header=True, index=False)
        os.system(f'chmod 775 {full_path}/pointfinder_summary_{report_time}.csv')
    elif type == 'amrfinder_mut':
        current_df = pd.concat([current_combined_df,new_df], sort=False)
        current_df.drop_duplicates(keep='last', inplace=True)
        os.system(f'mv {path_to_current_file} {full_path}/backup/')
        os.system(f'cp "{path_to_new_table}" {full_path}/backup/tables/{os.path.basename(path_to_new_table).replace(".csv",time.strftime("_%Y_%m_%d:%H:%M:%S.csv"))}')
        current_df.to_csv(f'{full_path}/amrfinder_mut_summary_{report_time}.csv', header=True, index=False)
        os.system(f'chmod 775 {full_path}/amrfinder_mut_summary_{report_time}.csv')


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
    res_df = process_report(os.path.abspath(parser.hamronization)) 
    try:
        rfp_df = pd.read_csv(os.path.abspath(parser.resfinder))
    except:
        rfp_df = pd.DataFrame(columns = [
            "sample_id","analysis_batch_id","# Antimicrobial","Class","WGS-predicted phenotype","Match","Genetic background"
        ])
    try:
        pf_df  = pd.read_csv(os.path.abspath(parser.pointfinder))
    except: 
        pf_df = pd.DataFrame(columns = [
            "sample_id","analysis_batch_id","Mutation","Nucleotide change","Amino acid change","Resistance","PMID"
        ])
    try:
        amfm_df  = pd.read_csv(os.path.abspath(parser.amrfinder_mut))
    except: 
        amfm_df = pd.DataFrame(columns = [
            "sample_id","analysis_batch_id",
            "Protein identifier","Contig id",
            "Start","Stop","Strand","Gene symbol",
            "Sequence name","Scope","Element type",
            "Element subtype","Class","Subclass","Method",
            "Target length","Reference sequence length",
            "% Coverage of reference sequence",
            "% Identity to reference sequence",
            "Alignment length","Accession of closest sequence",
            "Name of closest sequence","HMM id","HMM description"
        ])
    #update hamronization-aggregated summary table
    update_combined_resistance_file( 
        path_to_new_table=parser.hamronization ,
        current_combined_df=current_combined_tuple[1], 
        path_to_current_file=current_combined_tuple[0],
        new_df=res_df
    )

    #update resfinder-phenotype summary table
    update_combined_resistance_file(
        path_to_new_table=parser.resfinder, 
        path_to_current_file=current_combined_tuple[2], 
        current_combined_df=current_combined_tuple[3], 
        new_df=rfp_df,
        type="resfinder"
    )

    #update pointfinder summary table
    update_combined_resistance_file(
        path_to_new_table=parser.pointfinder, 
        path_to_current_file=current_combined_tuple[4], 
        current_combined_df=current_combined_tuple[5], 
        new_df=pf_df,
        type="pointfinder"
    )

    #update amrfinderplus mutation summary table
    update_combined_resistance_file(
        path_to_new_table=parser.amrfinder_mut, 
        path_to_current_file=current_combined_tuple[6], 
        current_combined_df=current_combined_tuple[7], 
        new_df=amfm_df,
        type="amrfinder_mut"
    )
