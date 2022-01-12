import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) #TO AVOID PANDAS CLUTTERING THE OUTPUT

import os, subprocess, re, time, sys, argparse, json, pandas as pd
from collections import deque
from pathlib import Path

full_path = Path(os.path.dirname(os.path.realpath(__file__))).parents[0] #KMERFINDER_SUITE DIR
parser = argparse.ArgumentParser(description='A script infer bacterial species from raw pair-end reads (kmerfinder).') #ARGPARSER OBJECT TO PROVIDE COMMAND-LINE FUNCTIONALITY
parser.add_argument('-t', '--top_hits', metavar='\b', help = 'Number of top kmerfinder hits to include in assembly report (default = 3).', default=3, required=False)
req_arg_grp = parser.add_argument_group('required arguments') #TO DISPLAY ARGUMENT UNDER REQUIRED HEADER IN HELP MESSAGE
req_arg_grp.add_argument('-n', '--num_files', metavar='\b', help = 'Number of samples to be processed in-parallel on different nodes of the cluster', default=None, required=True)

if len(sys.argv)==1: #IF NO COMMAND-LINE ARGUMENTS PROVIDED - DISPLAY HELP AND STOP SCRIPT EXCECUTION
    parser.print_help(sys.stderr)
    sys.exit(1)

#PARSE ARGUMENTS
args = parser.parse_args()
fastq_path = "/home/groups/nmrl/bact_analysis/kmerfinder_suite/input_files/" #PATH TO FASTQ FILES FROM COMMAND-LINE ARGUMENT
report_time = time.strftime("%m_%d_%Y") #DATESTAMP TO USE IN AGGREGATED REPORT
os.chdir("/home/groups/nmrl/bact_analysis/kmerfinder_suite/") #ALLOWS TO USE JOB SCRIPT FOR KMERFINDER

file_list = []
for (root,dirs,files) in os.walk(fastq_path, topdown=True): #COLLECTS FULL PATHS TO ALL FASTQ UNDER FASTQ_PATH AND SUBDIRS 
    for name in files: 
        if 'fastq.gz' in name: file_list.append(os.path.join(root, name)) 

if len(file_list) == 0: #NO FILES FOUND - NO PROCESSING NEEDED
    print('No fastq files were found in /home/groups/nmrl/bact_analysis/kmerfinder_suite/input_files/')
    sys.exit(1)


seq_id_list = set(file.split("/")[-1].split("_",1)[0] for file in file_list) #UNIQUE IDS EXTRACTED FROM EACH FILE
read_dict = {} #TO STORE FILE PATHS MAPPED TO ID
for id in seq_id_list: read_dict[id]=[] #MAP ID TO EMPTY LIST - ALTERNATIVELY - CHECK IF READ_DICT ALREADY HAS KEY - 2*LINEAR RUNTIME IN BOTH CASES, BUT THIS WAY CODE IS CLEANER
#MAP PATHS TO IDS
for id in seq_id_list: 
    for file_path in file_list:
        if id in file_path: #ID AS SUBSTRING IN PATH
            read_dict[id].append(file_path) #ADD FULL PATH TO FASTQ FILE TO DICT
            if len(read_dict[id]) == 2: #PATHS TO BOTH READ 1 & READ 2 FILES ARE IN THE DICT
                break #PROCEED TO THE NEXT ID

def submit_job(sample_id,que):
    '''Function to submit kmerfinder job on HPC'''
    read_dict[sample_id] = sorted(read_dict[sample_id]) #TO ENSURE THAT READ_1 IS SUPPLIED TO THE JOB SCRIPT BEFORE READ_2
    read_1 = read_dict[sample_id][0]
    read_2 = read_dict[sample_id][1]
    #ADDING VERBOUS MESSAGES REGARDING SUBMITTED JOBS
    print(sample_id, read_1, read_2)
    print(sample_id, f'{len(que)}/{len(list(read_dict.keys()))}')
    print(['qsub', '-F', f'{read_1} {read_2} {sample_id}', 'run_kmerfinder.sh'])
    #SYNTAX FOR QSUB WITH ARGS IS ... INTERESTING
    # subprocess.check_call(['qsub', '-F', f'{read_1} {read_2} {sample_id}', 'kmerfinder_job.sh'])

sample_limit = int(args.num_files) #NUMBER OF SAMPLES TO BE PROCESSED AT THE SAME TIME ON DIFFERENT NODES
sleep_time = 30 #IN SECONDS TO WAIT BEFORE ADDING NEW SAMPLES TO QUE - TO AVOID HAVING TOO MANY JOBS AT ONCE
sample_que = deque(read_dict.keys()) #CREATE A QUE OF SAMPLES WITHIN THE SCRIPT
#SUBMITTING JOBS
processing_over=False #CONTROL TO BREAK SUBMISSION IF JOB STACK IS FULL OR ALL SAMPLES ALREADY PROCESSED 
while len(sample_que) > 0: #UNTIL THERE ARE SAMPLES TO BE PROCESSED
    for _ in range(sample_limit):
        if len(sample_que) > 0:
            sample_id = sample_que.pop() #POP RIGHT-MOST SAMPLE ID FROM QUE (QUE IS SHORTENED BY ON ELEMENT)
            submit_job(sample_id,sample_que)
        else:
            processing_over=True
            break
    if processing_over:
        break
    time.sleep(sleep_time) #DELAY SUBMISSION 

file_dict = dict() #CHECK IF THE PROCESSING IS FINISHED BEFORE AGGREGATING - LOOKS FOR ALL JSON FILES UNDER KMERFINDER_SUITE/OUTPUT & SUBDIRS
while len(file_dict.keys()) < len(read_dict.keys()):
    for (root,dirs,files) in os.walk("/home/groups/nmrl/bact_analysis/kmerfinder_suite/output/", topdown=True): 
        for name in files: 
            jfile_path = os.path.join(root, name)
            if 'json' in jfile_path: 
                if os.path.dirname(jfile_path).split('/')[-1] in read_dict.keys(): #IF DIRECTORY NAME CORRESPONDS TO ID OF THE PROCESSED SAMPLE
                    file_dict[os.path.dirname(jfile_path).split('/')[-1]] = jfile_path #MAP PATH TO JSON FILE TO THE SAMPLE ID
    if len(file_dict.keys()) == len(read_dict.keys()): break #STOP WAITING IF ALL FILES WERE GENERATED
    print(f'Waiting for the processing to compete for {sleep_time} seconds.')
    time.sleep(sleep_time) #DELAY AGGREGATION IF SOME FILES ARE MISSING
os.system('rm find_bact.*') #REMOVING JOB REPORTS

#AGGREGATING RESULTS IN ONE TABLE
output_df = pd.DataFrame()
for id in read_dict.keys(): #ONLY FRESHLY PROCESSED SAMPLES
    with open(file_dict[id]) as json_file:
        data = json.load(json_file) #OPEN KMERFINDER REPORT
    hits = data['kmerfinder']['results']['species_hits'] #GET HITS DATA
    # CONVERT TO DATAFRAME
    sample_df = pd.DataFrame()
    for key in hits.keys(): #COMBINE HITS INTO DATAFRAME & ADD ID + ACCESSION
        row = hits[key]
        row['sample_id'] = id
        row['accession'] = key.split(" ",1)[0].split('.',1)[0]
        for key in row.keys(): #CONVERTING DATA TO PANDAS-ACCEPTED FORMAT
            row[key] = [row[key]]
        df = pd.DataFrame.from_dict(row)
        sample_df = sample_df.append(df)
    sample_df['Score'] = sample_df['Score'].astype(int)
    output_df = output_df.append(sample_df.sort_values('Score', ascending=False).head(int(args.top_hits))) #INCLUDE USER-DEFINED NUMBER OF HITS IN THE REPORT (OR ALL HITS, IF LESS THAN SPECIFIED IS AVAILABLE)

output_df.to_csv(f"/home/groups/nmrl/bact_analysis/kmerfinder_suite/output/{report_time}_assambled_kmerfinder_report.csv",header=True, index=False)