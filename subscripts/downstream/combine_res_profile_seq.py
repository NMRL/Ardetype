"""
2024-07-17
Script to run pairwise comparison of resistance markers identified by AMRFinder+, ResFinder and RGI.
"""

#########
# Imports
#########

import pandas as pd
import subprocess as sp
import os
import time
import re
import logging
import time
import functools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

########
# Setups
########

runtime_data = pd.DataFrame(columns=['function', 'runtime_sec', 'runtime_min'])
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
mash_container_path = "/mnt/beegfs2/home/groups/nmrl/image_files/mashtree_latest.sif"
bindash_container_path = "/mnt/beegfs2/home/groups/nmrl/bact_analysis/Ardetype/res_harm_seq/bindash_2.2--h43eeafb_0.sif"
dashing2_binary_path = "/mnt/beegfs2/home/groups/nmrl/bact_analysis/Ardetype/res_harm_seq/dashing2_savx2"
blast_container_path = "/mnt/beegfs2/home/groups/nmrl/bact_analysis/Ardetype/res_harm_seq/blast_latest.sif"

#####################
# Aggregation settings
#####################

tool_marks = ['RES', 'AMR+', 'RGI']
sample_identifier = 'sample_id'
specials = {
        "Escherichia_coli_AcrAB-TolC_with_MarR_mutations_conferring_resistance_to_ciprofloxacin_and_tetracycline": "AcrAB-TolC_with_MarR",
        "blaTEMp_G162T_Escherichia_amoxicillin-clavulanic_acid_piperacillin-tazobactam_ticarcillin-clavulanic_acid_resistant_blaTEM":"blaTEMp_G162T"
    }
common_columns = [
    "Reference-Tool",
    "Query-Tool",
    "Reference-Marker",
    "Query-Marker",
    'Contig-reference',
    'Contig-query',
    "Label_match"
]

mash_columns = [
    "Reference-ID",
    "Query-ID", 
    "Mash-distance",
    "P-value",
    "Matching-hashes",
]
mash_suffix = '_rmc_mash.csv'
mash_id = 'Mash-distance'
mash_threshold = 1
mash_sim = False

bindash_columns = [
    "Reference-ID",
    "Query-ID", 
    "BinDash-distance",
    "P-value",
    "Matching-hashes",
]
bindash_suffix = '_rmc_bindash.csv'
bindash_id = 'BinDash-distance'
bindash_threshold = 1
bindash_sim = False

dashing2_columns = [
    "Source",
    "Query", 
    "Similarity"
]
dashing2_suffix = '_rmc_dashing2.csv'
dashing2_id = 'Similarity'
dashing2_threshold = 0


blast_columns = [
        "query acc.ver", 
        "subject acc.ver", 
        "% identity", 
        "alignment length", 
        "mismatches", 
        "gap opens", 
        "q. start", 
        "q. end", 
        "s. start", 
        "s. end", 
        "evalue", 
        "bit score"
    ]
blast_suffix = '_rmc_blast.csv'
blast_id = '% identity'
blast_threshold = 0


##################
# Helper functions
##################

# Decorators

def measure_runtime(func):
    '''Decorator function to measure runtime of other functions in seconds
    and show it using logging module.
    '''
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        runtime = end_time - start_time
        if func.__name__ == "calculate_distances":
            name = kwargs['command_wrapper'].__name__
        else:
            name = func.__name__
        logger.info(f"Function '{name}' executed in {runtime:.4f} seconds")
        runtime_df = pd.DataFrame.from_dict(
            {
            'function'   : [name], 
            'runtime_sec': [runtime],
            'runtime_min': [runtime/60],
            }
        )
        global runtime_data
        runtime_data = pd.concat([runtime_data, runtime_df])
        return result
    return wrapper


def retry_on_failure(retries: int = 3, retry_delay: int = 10):
    '''
    Decorator to add retry logic to functions that may fail due to transient issues.
    
    Parameters:
        retries (int): Number of times to retry the function.
        retry_delay (int): Delay between retries in seconds.
    '''
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            for attempt in range(retries):
                # Call the wrapped function and capture its result
                result = func(*args, **kwargs)
                # Check for errors in stderr
                error_signatures = ["fatal", "error", "segmentation"]
                stderr = result[2].lower()
                return_code = int(result[-1])
                error_signs = any([es in stderr for es in error_signatures])
                if not error_signs and return_code == 0:
                    return result
                else:
                    if attempt < retries - 1:
                        time.sleep(retry_delay)
            # Return the last result after all retries
            return result
        return wrapper
    return decorator

# Plain functions

def write_fasta(records:list, path:str):
    '''
    Oneliner:
        write list of biopython SeqRecord objects
        into fasta file using SeqIO wrapper.
    '''
    with open(path, 'w+') as f:
        SeqIO.write(records, f, "fasta")


@measure_runtime
def combine_profiles(path_list:list, out_path:str):
    '''
    Oneliner:
        sequences aggregation and writing into single fasta file.
    '''
    combined_records = combine_fasta(path_list)
    write_fasta(combined_records, out_path)


@measure_runtime
def generate_split(file_path:str, split_folder_path:str):
    '''
    Bash helper to create a folder <split_folder_path> 
    containing all sequences from fasta <file_path> as individual fasta files.
    '''
    try:
        os.makedirs(split_folder_path, exist_ok=False)
        bash_script = f'''
        output_dir={split_folder_path}
        fasta={file_path}
        # Use awk to split multifasta into individual files with header names
        awk '/^>/{{header=$0; gsub(/[[:space:]]/, "_", header); sub(">", "", header); outname=sprintf("%s/%s.fasta", "'$output_dir'", header); print > outname; next;}} {{if(length($0)>0)print >> outname;}}' "$fasta"

        echo "Multifasta file '{file_path}' split into individual sequence files in '{split_folder_path}' - $(ls {split_folder_path} | grep .fasta | wc -l) sequences in total."
            '''
        output = sp.run(bash_script, shell=True, capture_output=True)
        logger.info(output.stdout.decode('utf-8'))
    except FileExistsError:
        logger.info(f"Skipping fasta split generation - folder already exists.")


def get_sketch_combinations(
    split_folder_path:str, 
    sketch_tag:str, 
    tool_marks:list) -> list:
    '''
    Wrapper:
        Given file extension to detect single tool-specific sketches as <sketch_tag>
        and tool identifiers (see combine_fasta function) as <tool_marks>,
        generates a list of unique combinations of sequences, 
        excluding self-with-self combinations based on tool marks.
    '''
    sketch_list = [os.path.join(split_folder_path,p) for p in os.listdir(split_folder_path) if p.endswith(sketch_tag)]
    if len(sketch_list) > 0:
        combs = list(combinations(sketch_list, 2))
        filtered_combs = [c for c in combs if not any(mark in c[0] and mark in c[1] for mark in tool_marks)]
        return sorted(filtered_combs)
    else:
        return []


def parallel_sketch_comparison(
    command_wrapper, 
    filtered_combs:list, 
    split_folder_path:str,
    dist_file_path:str,
    container_path:str,
    max_workers:int=1):
    '''
    Wrapper:
        Given function that wrapps bash command for distance calculation as <command_wrapper>*,
        list of sketch combinations to check as <filtered_combs> (see get_sketch_combinations),
        path to the folder where resistance markers are saved in separate fasta files as <split_folder_path>(see generate_split)
        path to file where pairwise distances should be aggregated as <dist_file_path>,
        and path to the container for the tool used for pairwise distance calculation as <container_path,
        runs the <command_wrapper> on multiple processors using concurrent.futures.ProcessPoolExecutor utility.
        Returns list of dict containing information about executed commands, stdout and stderr (convertable to pandas DataFrame).
        *expects to return the tuple of 3 strings - command, stdout and stderr (see execute_mash_command for example).
    '''
    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(command_wrapper, c, split_folder_path, dist_file_path, container_path) for c in filtered_combs]
        for future in futures:
            command, stdout, stderr, return_code = future.result()
            results.append({'command': command, 'stdout': stdout, 'stderr': stderr, 'return_code':return_code})
    return results


def save_comparison_logs(logs:list, log_path:str):
    '''
    Given list of dict containing command, stdout, stderr string triplets as <logs> (see <command_wrapper> examples)
    converts it to dataframe, aggregates by stderr to see distinct error counts,
    and saves as csv file at <log_path>.
    '''
    logs = pd.DataFrame(logs)
    logs = logs.groupby(['stderr', 'return_code']).count().reset_index()[['stderr', 'return_code', 'command']]
    logs = logs.rename(columns={'command':'entry_count'})
    logs.to_csv(log_path, header=True, index=False)


def normalize_columns(df:pd.DataFrame, reference:str, query:str):
    '''
    Given dataframe of pairwise comparisons as <df>,
    name of the reference sequence column as <reference> and query sequence columns as <query>,
    generates 6 tool-agnostic columns based on fasta header structure introduced by combine_fasta function.
    '''
    df[reference]          = df[reference].apply(lambda x: os.path.basename(x))
    df[query]              = df[query].apply(lambda x: os.path.basename(x))
    df["Reference-Tool"]   = df[reference].str.split('|', expand=True)[2].str.replace(r'.fasta', '')
    df["Query-Tool"]       = df[query].str.split('|', expand=True)[2].str.replace(r'.fasta', '')
    df["Reference-Marker"] = df[reference].str.split('|', n=2, expand=True)[0]
    df["Query-Marker"]     = df[query].str.split('|', n=2, expand=True)[0]
    df["Contig-reference"] = df[reference].str.split('|', n=2, expand=True)[1]
    df["Contig-query"]     = df[query].str.split('|', n=2, expand=True)[1]
    return df


def match_marker_labels(df:pd.DataFrame, qry_mark:str, ref_mark:str):
    '''
    Given dataframe <df> containing normalized resistance marker info as <qry_mark> and <ref_mark> columns,
    constructs attempts to check if both columns contain the same marker names, using simple splitting logic.
    Returns list of booleans (0/1) to be appended as column to <df> (matched in length to its index).
    '''
    tp = []
    for i, row in df.iterrows():
        # Try extracting marker names
        norm_q_mark = row[qry_mark].lower()
        norm_q_mark = norm_q_mark.split('_',maxsplit=1)[0]
        norm_r_mark = row[ref_mark].lower()
        norm_r_mark = norm_r_mark.split('_',maxsplit=1)[0]
        marker_name = min(norm_q_mark, norm_r_mark, key=len)
        
        # Match conditions based on name presence
        in_query = marker_name in row[qry_mark].lower()
        in_ref   = marker_name in row[ref_mark].lower()
        match_condition = in_query and in_ref
        if match_condition:
            tp.append(1)
        else:
            tp.append(0)
    return tp


@measure_runtime
def calculate_distances(
    command_wrapper,
    split_folder_path:str, 
    dist_file_path:str, 
    container_path:str,
    tool_marks:list,
    log_path:str,
    max_workers:int=10,
    sketch_tag:str='.msh',
    dist_file_sep:str='\t',
    tool_columns:list=mash_columns,
    normalized_columns:list=common_columns,
    reference_column="Reference-ID",
    query_column="Query-ID"):
    '''
    Aggregator:
        Combination list construction and parallel (not working) aggregation 
        of pairwise Jackard distances based on set of sketches.
    
        Logging to be implemented based on dashing2 template.
    '''
    filtered_combs = get_sketch_combinations(split_folder_path, sketch_tag=sketch_tag, tool_marks=tool_marks)
    if len(filtered_combs) > 0:
        if os.path.isfile(dist_file_path): os.remove(dist_file_path)
        logs = parallel_sketch_comparison(
            command_wrapper   = command_wrapper, 
            filtered_combs    = filtered_combs,
            split_folder_path = split_folder_path,
            dist_file_path    = dist_file_path,
            container_path    = container_path,
            max_workers       = max_workers
            )
        
        save_comparison_logs(logs, log_path)

        if os.path.isfile(dist_file_path) and os.stat(dist_file_path).st_size != 0:
            df = pd.read_csv(dist_file_path, sep=dist_file_sep, header=None)
            if len(df) > 0:
                df.columns = tool_columns
                df = df.drop_duplicates([reference_column, query_column]) # parallel processing retries may cause duplications
                df = normalize_columns(df, reference_column, query_column)
                label_matches = match_marker_labels(df, 
                    ref_mark = normalized_columns[2], 
                    qry_mark = normalized_columns[3])
                df[normalized_columns[-1]] = label_matches
            else:
                df = pd.DataFrame(columns = tool_columns + normalized_columns)
        else:
            df = pd.DataFrame(columns = tool_columns + normalized_columns)
    else:
        df = pd.DataFrame(columns = tool_columns + normalized_columns)
    df.to_csv(dist_file_path, header=True, index=False)


####################
# Specific functions
####################

@measure_runtime
def rgi_to_fasta(
    path:str, 
    head_columns:list=['Contig','Start', 'Stop', 'Best_Hit_ARO', 'Predicted_DNA'], 
    sep:str='\t'):
    '''
    Format changer: rgi.txt > fasta 
    Fasta structure in terms of biopython SeqRecord attributes:
        id: 'Contig':'Start'-'Stop'
        description: 'Best_Hit_ARO'
        seq: 'Predicted_DNA'
    '''

    df = pd.read_csv(path, sep=sep)
    df = df[head_columns]
    df['id'] = df['Contig'] + ':' + df['Start'].astype(str) + '-' + df['Stop'].astype(str)
    df['id'] = df['id'].str.replace(' ', '_')
    seqio_records = [SeqRecord(Seq(row['Predicted_DNA']), id=row['id'], description=row["Best_Hit_ARO"]) for _, row in df.iterrows()]
    return seqio_records


def read_fasta(path:str):
    '''Oneliner:
        Reads fasta to list of SeqRecord objects using biopython SeqIO parser.'''
    with open(path, 'r') as f:
        records = list(SeqIO.parse(f, 'fasta'))
    return records


@measure_runtime
def combine_fasta(
    path_list:list, 
    specials:dict=specials):
    '''
    Aggregator:
        reading, header format unification and
        aggregation of data from fasta files of <path_list> generated as
        output of resistance marker detections tools
        into single list of biopython SeqRecord objects.
        <specials> define special string compression rules to be applied to avoid 
        Exceptions due too long file names.
    '''
    combined_records = [] #headers + sequences
    current_tool     = '' #name of the tool for logging
    counts           = [] #number of records per tool

    #Header processing
    for path in path_list:
        if 'rgi' in path:
            current_tool = 'RGI'
            records      = rgi_to_fasta(path)
            for record in records:
                record.id = f"{record.description.replace('/','_')}|{record.id.split(':')[0]}"
                if record.id.endswith('_'):
                    record.id = record.id[:-1]
                record.id         += '|RGI'
                record.id          = record.id.replace(' ', '_') 
                record.description = ""
                for s in specials:
                    if s in record.id:
                        record.id = record.id.replace(s, specials[s])
        else:
            records = read_fasta(path)
            if 'ResFinder' in path:
                current_tool = 'RES'
                for record in records:
                    try:
                        record.description = record.description.replace('/','_')
                        record.description = record.description.split(', ')[5]
                        record.description = record.description.split(' len')[0]
                    except IndexError:
                        record.description = record.description.replace(' ', '_')
                    record.id = f"{record.id.replace(' ','_')}{record.description}|RES"
                    record.id = record.id.replace(',Contig name: ', '|')
                    record.description = ''
                    for s in specials:
                        if s in record.id:
                            record.id = record.id.replace(s, specials[s])
            elif 'amrfinderplus' in path:
                current_tool = 'AMR+'
                for record in records:
                    record.description = record.description.replace("/","_")
                    record.description = record.description.split(' ', maxsplit=2)[-1]
                    record.id = f"{record.description}|{record.id.split(':')[0]}"
                    record.id         += '|AMR+'
                    record.id = record.id.replace(' ', '_')
                    record.description = ""
                    for s in specials:
                        if s in record.id:
                            record.id = record.id.replace(s, specials[s])
        #Aggregation and counting
        seq_count = len(records)
        unq_seq_count = len(set(r.id for r in records))
        counts.append(unq_seq_count)
        combined_records.extend(records)
        logger.info(f'Tool: {current_tool} - {seq_count} sequences found, {unq_seq_count} unique headers.')

    #Final aggregation and numbers
    total_records = len(combined_records)
    comb_estimates = counts[0]*counts[1] + counts[0]*counts[2] + counts[1]*counts[2]
    logger.info(f'Total number of sequences found - {total_records}.')
    logger.info(f'Estimated number of combinations - {comb_estimates}.')
    return combined_records


##################
# Command wrappers
##################


@measure_runtime
def generate_mash_sketches(
    split_folder_path:str, 
    container_path: str, 
    log_path:str='mash_sketch.log'):
    '''
    Bash helper to generate mash sketches from fasta files in <split_folder_path>;
    '''

    bash_script = f'''
    module load singularity
    for f in {split_folder_path}/*.fasta; do
        singularity run --bind {split_folder_path}:{split_folder_path} \
        {container_path} mash sketch -o "${{f/.fasta/.msh}}" "${{f}}"
    done
    '''
    output = sp.run(bash_script, shell=True, capture_output=True)
    logging.info(f'Finished mash sketch generation for - {split_folder_path}.')
    with open(log_path, 'w') as f:
        f.write(output.stdout.decode("utf-8"))
        f.write(output.stderr.decode("utf-8"))


@measure_runtime
def generate_bindash_sketches(
    split_folder_path:str, 
    container_path:str, 
    log_path:str='bindash_sketch.log'):
    '''
    Bash helper to generate bindash sketches from fasta files in <split_folder_path>;
    '''
    bash_script = f'''
    module load singularity
    for f in {split_folder_path}/*.fasta;  do 
        singularity run --bind {split_folder_path}:{split_folder_path} \
        {container_path} bindash sketch --minhashtype=-1 --outfname=${{f/.fasta/.bsh}} ${{f}}; 
    done
    '''
    output = sp.run(bash_script, shell=True, capture_output=True)
    logging.info(f'Finished bindash sketch generation for - {split_folder_path}.')
    with open(log_path, 'w') as f:
        f.write(output.stdout.decode("utf-8"))
        f.write(output.stderr.decode("utf-8"))


@retry_on_failure(retries=25, retry_delay=10)
def execute_mash_command(
    c: tuple, 
    split_folder_path: str, 
    dist_file_path: str,
    container_path: str
) -> tuple:
    '''
    Bash helper: mash distance calculation command;
    '''
    command = f'''
    module load singularity
    singularity run --bind \
        {split_folder_path}:{split_folder_path} \
        {container_path} mash dist \
        "{c[0]}" "{c[1]}" >> {dist_file_path}
    '''
    result = sp.run(command, shell=True, capture_output=True)
    stderr = result.stderr.decode('utf-8')
    stdout = result.stdout.decode('utf-8')
    return command, stdout, stderr, result.returncode


@retry_on_failure(retries=25, retry_delay=10)
def execute_bindash_command(
    c: tuple, 
    split_folder_path: str, 
    dist_file_path: str,
    container_path: str
) -> tuple:
    '''
    Bash helper: bindash distance calculation command;
    Retries with delay added to deal with race conditions;
    '''
    # Construct the output file name
    c1 = c[0].replace('.bsh','')
    c2 = os.path.basename(c[1]).replace('.bsh','')
    outfname = f"{c1}_{c2}.bindash"
    # Construct the command
    command = f'''
    module load singularity
    singularity run --bind {split_folder_path}:{split_folder_path} \
        {container_path} bindash dist --outfname="{outfname}"\
        "{c[0]}" "{c[1]}"
    cat "{outfname}" >> "{dist_file_path}"
    '''
    result = sp.run(command, shell=True, capture_output=True)
    stderr = result.stderr.decode('utf-8')
    stdout = result.stdout.decode('utf-8')
    return command, stdout, stderr, result.returncode


@retry_on_failure(retries=3, retry_delay=2)
def execute_dashing2_command(
    c: tuple, 
    split_folder_path: str, 
    dist_file_path: str,
    container_path: str
) -> tuple:
    '''
    Bash helper: dashing2 distance calculation command;
    '''

    # Construct the output file name
    c1 = c[0].replace('.fasta','')
    c2 = os.path.basename(c[1]).replace('.fasta','')
    outfname = f"{c1}_{c2}.d2"
    command = f'{container_path} sketch --cmpout "{outfname}" "{c[0]}" "{c[1]}"'
    result = sp.run(command, shell=True, capture_output=True)
    stderr = result.stderr.decode('utf-8')
    stdout = result.stdout.decode('utf-8')
    
    with open(outfname, 'r') as f:
        contents = f.readlines()
    contents = contents[-3:]
    for i,s in enumerate(contents):
        contents[i] = s.strip().split()
        if contents[i].count('-') == 2:
            contents.pop(i)
    empty_index = contents[1].index('-')
    for i,s in enumerate(contents):
        contents[i].pop(empty_index)
    data_string = f'"{contents[0][1]}","{contents[1][0]}","{contents[1][1]}"\n'
    with open(outfname, 'w') as f:
        f.write(data_string)

    os.system(f'cat "{outfname}" >> {dist_file_path}')
    return command, stdout, stderr, result.returncode


@retry_on_failure(retries=25, retry_delay=45)
def execute_blast_command(
    c: tuple, 
    split_folder_path: str, 
    dist_file_path: str,
    container_path: str
) -> tuple:
    '''
    Bash helper: blast distance calculation command;
    '''
    # Construct the output file name
    c1 = c[0].replace('.fasta','')
    c2 = os.path.basename(c[1]).replace('.fasta','')
    outfname = f"{c1}_{c2}.blast"
    command = f'''
    module load singularity
    singularity run --bind \
        {split_folder_path}:{split_folder_path} \
            {container_path} blastn -query "{c[0]}" -subject "{c[1]}" -out "{outfname}" -outfmt 7'''
    result = sp.run(command, shell=True, capture_output=True)
   
    if os.path.isfile(outfname):
        with open(outfname, 'r') as f:
            contents = f.readlines()
        try:
            if not contents:
                data = [0 for _ in blast_columns]
            elif '0 hits' in contents[3]:
                data = [0 for _ in blast_columns]
            else:
                data = [v.strip() for v in contents[5].split('\t')]
            data[0] = os.path.basename(c1)
            data[1] = c2
            data_string = ",".join([f'"{e}"' for e in data])
            data_string += "\n"
            with open(outfname, 'w') as f:
                f.write(data_string)
        except IndexError as e:
            stderr = f'ERROR attempting to parse {outfname}'
            stdout = ""
            return_code = -1
            return command, stdout, stderr, return_code

        os.system(f'cat "{outfname}" >> {dist_file_path}')

    stderr = result.stderr.decode('utf-8')
    stdout = result.stdout.decode('utf-8')
    return command, stdout, stderr, result.returncode


@measure_runtime
def main(
    path_list, #RGI, RF, AMR+
    out_path, #combined fasta path
    split_folder_path, #Split fasta folder path
    mash_df_path,
    mash_sketch_log,
    mash_dist_log,
    dash2_df_path, 
    dash2_log_path,
    blast_df_path,
    blast_log_path,
    bindash_df_path, 
    bindash_sketch_log,
    bindash_distance_log,
    runtime_est_path
    ):
    combine_profiles(path_list, out_path)
    generate_split(out_path, split_folder_path)
    generate_bindash_sketches(split_folder_path, container_path=bindash_container_path, log_path=bindash_sketch_log)
    calculate_distances(
            command_wrapper=execute_bindash_command,
            split_folder_path=split_folder_path, 
            dist_file_path=bindash_df_path, 
            container_path=bindash_container_path,
            tool_marks=tool_marks,
            log_path=bindash_distance_log,
            max_workers=30,
            sketch_tag='.bsh',
            dist_file_sep='\t',
            tool_columns=bindash_columns,
            normalized_columns=common_columns,
            reference_column=bindash_columns[0],
            query_column=bindash_columns[1])
    generate_mash_sketches(
        split_folder_path=split_folder_path, 
        container_path=mash_container_path, 
        log_path=mash_sketch_log)
    calculate_distances(
        command_wrapper=execute_mash_command,
        split_folder_path=split_folder_path, 
        dist_file_path=mash_df_path, 
        container_path=mash_container_path,
        tool_marks=tool_marks,
        log_path=mash_dist_log,
        max_workers=10,
        sketch_tag='.msh',
        dist_file_sep='\t',
        tool_columns=mash_columns,
        normalized_columns=common_columns,
        reference_column=mash_columns[0],
        query_column=mash_columns[1])

    calculate_distances(
        command_wrapper=execute_dashing2_command,
        split_folder_path=split_folder_path, 
        dist_file_path=dash2_df_path, 
        container_path=dashing2_binary_path,
        tool_marks=tool_marks,
        log_path=dash2_log_path,
        max_workers=30,
        sketch_tag='.fasta',
        dist_file_sep=',',
        tool_columns=dashing2_columns,
        normalized_columns=common_columns,
        reference_column=dashing2_columns[0],
        query_column=dashing2_columns[1])

    calculate_distances(
        command_wrapper=execute_blast_command,
        split_folder_path=split_folder_path, 
        dist_file_path=blast_df_path, 
        container_path=blast_container_path,
        tool_marks=tool_marks,
        log_path=blast_log_path,
        max_workers=30,
        sketch_tag='.fasta',
        dist_file_sep=',',
        tool_columns=blast_columns,
        normalized_columns=common_columns,
        reference_column=blast_columns[0],
        query_column=blast_columns[1])
    runtime_data.to_csv(runtime_est_path, header=True, index=False)


##################
# Aggregator logic
##################

def batch_runtime(path_list):
    '''Aggregator:
        Calculates average & SD of runtime per number of res.marker combinations from files in <path_list>.
        
    long code lines
    review correctnes of calculations
    depends of file formats to work without exposing them as parameters
    '''
    from numpy import inf
    runtime_aggr = pd.DataFrame(columns = ["function","runtime_sec","runtime_min","path", "comb_count", "norm_runtime_sec"])
    for path in path_list:
        df = pd.read_csv(path)
        df['path'] = [path for _ in df.index]
        blast = pd.read_csv(path.replace('_rmc_runtime.csv','_rmc_blast.csv'))
        df['comb_count'] = [len(blast) for _ in df.index]
        df['norm_runtime_sec'] = df['runtime_sec']/df['comb_count']
        df['norm_runtime_sec'] = df['norm_runtime_sec'].replace(inf, 0)
        runtime_aggr = pd.concat([runtime_aggr, df])
    runtime_aggr = runtime_aggr.groupby(['function']).agg({'norm_runtime_sec':('mean', 'std')})['norm_runtime_sec'].reset_index().rename(columns={'mean':'relative_runtime_avg_sec', 'std':'relative_runtime_std_sec'}).sort_values(by='relative_runtime_avg_sec')
    return runtime_aggr


def batch_tool(
    path_list, 
    columns, 
    suffix, 
    id_column, 
    id_threshold, 
    sample_id, 
    is_similarity=True):
    '''
    Aggregator:
        aggregates pairs identified by specific tool for all samples in the batch.
    '''
    full_columns = [sample_id]+columns 
    gt = pd.DataFrame(columns=columns)
    for path in path_list:
        df = pd.read_csv(path)
        if len(df) > 0:
            if is_similarity:
                df = df[df[id_column] > id_threshold]
            else:
                df = df[df[id_column] < id_threshold]
            df[sample_id] = [os.path.basename(path).replace(suffix,'') for _ in df.index]
            df = df[full_columns]
            gt = pd.concat([gt,df])
    gt = gt.sort_values(by=sample_id)
    return gt


def batch_test(path_list:list, gt:pd.DataFrame, gt_control_col:str, gt_control_metric:int, right_keys:list, left_keys:list, control_columns:list, dist_sim_def:list):
    '''Aggregator:
        performs merge-based comparison of pairs detected by tools (as files specified in <path_list>)
        and ground truth, supplied as dataframe (gt)
    
    '''
    gt_match_fraction = pd.DataFrame(columns = [
        'report',
        'false_negative',
        'false_positives'
    ])
    test_details = []
    for path, id_col, left_key, right_key, metric in zip(path_list, control_columns, left_keys, right_keys, dist_sim_def):
        df = pd.read_csv(path)
        for k in right_key: df[k] = df[k].str.replace('.fasta', '')
        if gt_control_metric > 0: gtf = gt[gt[gt_control_col] < gt_control_metric]
        else: gtf = gt[gt[gt_control_col] > gt_control_metric]
        if metric > 0: dff = df[df[id_col] < metric]
        else: dff = df[df[id_col] > metric]

        #Get sets of all label pair combinations for ground truth & query (both ways, because order of comparison is not guaranteed but distance is symmetric)
        gt_set = set((gtf[left_key[0]]+"||"+gtf[left_key[1]]).to_list() + (gtf[left_key[1]]+"||"+gtf[left_key[0]]).to_list())
        check_set = set((dff[right_key[0]]+"||"+dff[right_key[1]]).to_list() + (dff[right_key[1]]+"||"+dff[right_key[0]]).to_list())

        fn = [[],[]]
        #Get labels of false negatives as labels unique to ground truth
        fn_set = gt_set.difference(gt_set.intersection(check_set))
        for key in fn_set:
            split = key.split("||")
            fn[0].append(split[0])
            fn[1].append(split[1])

        #Filter records corresponding to false negatives based on constructed key set
        fn_df = pd.DataFrame()
        if fn[0]:
            for e1,e2 in zip(fn[0], fn[1]):
                tdf = gt[(gt[left_key[0]] == e1) & (gt[left_key[1]] == e2)]
                fn_df = pd.concat([fn_df, tdf])
        fn_count = len(fn_df)

        # Get labels of false positives as labels unique to query
        fp = [[],[]]
        fp_set = check_set.difference(gt_set.intersection(check_set))
        for key in fp_set:
            split = key.split("||")
            fp[0].append(split[0])
            fp[1].append(split[1])

        #Filter records corresponding to false positives based on constructed key set
        fp_df = pd.DataFrame()
        if fp[0]:
            for e1,e2 in zip(fp[0], fp[1]):
                tdf = df[(df[right_key[0]] == e1) & (df[right_key[1]] == e2)]
                fp_df = pd.concat([fp_df, tdf])
        fp_count = len(fp_df)

        #Aggregate the records
        test_details.append(pd.concat([fn_df,fp_df]))

        #Aggregate the counts
        test_results = pd.DataFrame.from_dict({gt_match_fraction.columns[0]:[os.path.basename(path)], gt_match_fraction.columns[1] : [fn_count], gt_match_fraction.columns[2]:[fp_count]})
        gt_match_fraction = pd.concat([gt_match_fraction, test_results])
    return gt_match_fraction, test_details


#################
# Consensus logic
#################

def aggregate_distances(csv_path, sample_id):
    # Read the CSV file
    df = pd.read_csv(csv_path)
    
    # Define columns for aggregation
    agg_columns = ['sample_id','gene_symbol', 'input_sequence_id', 'db_vote', 'jackard_similarities', 'db_sources', 'original_gene_name']
    
    # Initialize results DataFrame
    results = pd.DataFrame()
    
    # Separate rows based on Similarity
    multiples = df[df['Similarity'] != 0]
    singles = df[df['Similarity'] == 0]
    unq_singles = set(singles['Reference-Marker'].to_list() + singles['Query-Marker'].to_list())
    
    # Process rows with Similarity == 0
    for s in unq_singles:
        tools_ref = singles[(singles['Reference-Marker'] == s)]['Reference-Tool']
        contigs_ref = singles[(singles['Reference-Marker'] == s)]['Contig-reference']
        tools_qry = singles[(singles['Query-Marker'] == s)]['Query-Tool']
        contigs_qry = singles[(singles['Query-Marker'] == s)]['Contig-query']
        tools = list(set(tools_ref.to_list() + tools_qry.to_list()))
        contigs = list(set(c.split("_")[0] for c in contigs_ref.to_list() + contigs_qry.to_list()))
        dff = pd.DataFrame(dict(zip(agg_columns, [[sample_id],[s.lower()], [contigs[0]], [1], ['NA'], [tools[0]], [s]])))
        results = pd.concat([results, dff])
    
    # Process rows with Similarity != 0
    unq_mult_markers = set(multiples['Reference-Marker'].to_list() + multiples['Query-Marker'].to_list())
    vote_counts = {}
    for m in unq_mult_markers:
        vote_counts_ref = multiples[(multiples['Reference-Marker'] == m)]
        vote_counts_qry = multiples[(multiples['Query-Marker'] == m)]
        dd = pd.concat([vote_counts_ref, vote_counts_qry])
        original_names = list(set(dd['Reference-Marker'].to_list() + dd['Query-Marker'].to_list()))
        
        contigs_r = multiples[(multiples['Reference-Marker'] == m)]['Contig-reference']
        contigs_rq = multiples[(multiples['Reference-Marker'] == m)]['Contig-query']
        contigs_q = multiples[(multiples['Query-Marker'] == m)]['Contig-query']
        contigs_qr = multiples[(multiples['Query-Marker'] == m)]['Contig-reference']
        contigs = list(set(c.split("_")[0] for c in contigs_r.to_list() + contigs_q.to_list() + contigs_rq.to_list() + contigs_qr.to_list()))
        
        tools_r = multiples[(multiples['Reference-Marker'] == m)]['Reference-Tool']
        tools_rq = multiples[(multiples['Reference-Marker'] == m)]['Query-Tool']
        tools_q = multiples[(multiples['Query-Marker'] == m)]['Query-Tool']
        tools_qr = multiples[(multiples['Query-Marker'] == m)]['Reference-Tool']
        
        simil_r = multiples[(multiples['Reference-Marker'] == m)]['Similarity']
        simil_q = multiples[(multiples['Query-Marker'] == m)]['Similarity']
        
        similarity = sorted([str(i) for i in list(set(simil_r.to_list() + simil_q.to_list()))])
        tools = list(set(tools_r.to_list() + tools_q.to_list() + tools_rq.to_list() + tools_qr.to_list()))
        vote_count = len(tools)
        dff = pd.DataFrame(dict(zip(agg_columns, [[sample_id], [m.lower()], [" | ".join(contigs)], [vote_count], [" | ".join(similarity)], [" | ".join(tools)], [" | ".join(original_names)]])))
        results = pd.concat([results, dff]).astype(str)
        
        df_sorted = results.sort_values(by=['gene_symbol', 'db_vote'], ascending=[True, True])
        results = df_sorted.drop_duplicates(subset='gene_symbol', keep='last')
        results = results.drop_duplicates(subset='original_gene_name', keep='first')
    return results


if __name__ == "__main__":
    path_list = [
        '/mnt/home/jevgen01/nmrl/bact_analysis/Ardetype/240614_NDX550703_RUO_0109_AHHTNWBGXW_20240620_162245/folded_2406013429_S84_output/2406013429_S84.rgi.txt',
        '/mnt/home/jevgen01/nmrl/bact_analysis/Ardetype/240614_NDX550703_RUO_0109_AHHTNWBGXW_20240620_162245/folded_2406013429_S84_output/2406013429_S84_resfinder/ResFinder_Hit_in_genome_seq.fsa', #ResFinder_Resistance_gene_seq.fsa',
        '/mnt/home/jevgen01/nmrl/bact_analysis/Ardetype/240614_NDX550703_RUO_0109_AHHTNWBGXW_20240620_162245/folded_2406013429_S84_output/2406013429_S84_amrfinderplus_found_seq.fasta'
        ]
    out_path = 'test.fasta'
    split_folder_path = '/mnt/beegfs2/home/groups/nmrl/bact_analysis/Ardetype/res_harm_seq/split_test'
    mash_df_path = 'distances_mash.csv'
    dash2_df_path = 'distances_dashing2.csv'
    blast_df_path = 'blast_results.csv'
    blast_log_path = 'blast_log.csv'
    dash2_log_path = 'dashing2_log.csv'
    mash_sketch_log = 'mash_sketch.log'
    mash_distance_log = 'mash_dist_log.csv'
    bindash_df_path = 'distances_binhash.csv'
    bindash_sketch_log = 'bindash_sketch.log'
    bindash_distance_log = 'bindash_dist_log.csv'
    runtime_est_path = 'runtime_estimates.csv'
    # main(
    #     path_list,
    #     out_path,
    #     split_folder_path,
    #     mash_df_path,
    #     mash_sketch_log,
    #     mash_distance_log,
    #     dash2_df_path,
    #     dash2_log_path,
    #     blast_df_path,
    #     blast_log_path,
    #     bindash_df_path,
    #     bindash_sketch_log,
    #     bindash_distance_log,
    #     runtime_est_path
    # )
    blast = pd.read_csv(blast_df_path)

    tool_test, merges = batch_test(
        path_list = [mash_df_path, bindash_df_path, dash2_df_path],
        gt = blast,
        gt_control_col = '% identity',
        gt_control_metric = 0,
        right_keys = [["Reference-ID", "Query-ID"]]*2 + [["Query", "Source"]],
        left_keys = [["query acc.ver","subject acc.ver"]]*3, 
        control_columns=[mash_id, bindash_id, dashing2_id],
        dist_sim_def=[1,1,0])
    for r, name in zip(merges, ['mash', 'bindash', 'dashing2']):
        r.to_csv(f'rmc_benchmark_details_{name}.csv', header=True, index=False)
    tool_test.to_csv("rmc_benchmark_report.csv", header=True, index=False)
    
    #Consensus
    input_file = '/mnt/home/jevgen01/nmrl/bact_analysis/Ardetype/res_harm_seq/distances_dashing2.csv'
    output_file = 'dashing2_consensus.csv'
    agg = aggregate_distances(input_file, sample_id="22364225_S21")
    agg.to_csv(output_file,index=False)
