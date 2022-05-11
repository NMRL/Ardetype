import os, sys, yaml, subprocess, pandas as pd, shutil, time, re, argparse
from getpass import getuser

def parse_folder(folder_pth_str, file_fmt_str, substr_lst=None, regstr_lst=None):
    '''
    Given path to the folder (folder_pth_str) and file format (file_fmt_str), returns a list, 
    containing absolute paths to all files of specified format found in folder and subfolders,
    except for files that contain patterns to exclude (specified in regstr_lst) or substrings to exclude (specified in substr_lst).    
    '''
    name_series = pd.Series(dtype="str") #initialize pandas series to store path values
    for (root,dirs,files) in os.walk(folder_pth_str, topdown=True): #get list of file paths (from parent dir & subdir)
        new_files = pd.Series(files, dtype="str") #convert file names in new folder to pandas series
        new_files = new_files[new_files.str.contains(file_fmt_str)] #keep only paths to files of specifed format
        new_files = f"{os.path.abspath(root)}/"  + new_files #append absolute path to the file
        if substr_lst is not None and regstr_lst is not None and (len(substr_lst) + len(substr_lst) > 0): #if both regex and substring filters provided
            new_files = new_files[~new_files.str.contains('|'.join(substr_lst+regstr_lst))].reset_index(drop=True)
        elif regstr_lst is not None and len(regstr_lst) > 0: #if only regex filters provided
            if len(regstr_lst) > 1: #checking single filter case
                new_files = new_files[~new_files.str.contains('|'.join(regstr_lst))].reset_index(drop=True)
            else:
                new_files = new_files[~new_files.str.contains(regstr_lst[0])].reset_index(drop=True)
        elif substr_lst is not None and len(substr_lst) > 0: #if only substring filters provided
            if len(substr_lst) > 1: #checking single filter case
                new_files = new_files[~new_files.str.contains('|'.join(substr_lst))].reset_index(drop=True)
            else:
                new_files = new_files[~new_files.str.contains(substr_lst[0])].reset_index(drop=True)
        name_series = name_series.append(new_files).reset_index(drop=True) #aggregating filtered paths
    return name_series.tolist()


def create_sample_sheet(file_lst, generic_str, regex_str=None, mode=0):
    """
    Given (list) of paths to files and a generic part of the file name (e.g. _contigs.fasta or _R[1,2]_001.fastq.gz string, regex expected for fastq), mode value (int 1 for fasta, 0 (default) for fastq)
    and a sample_id regex pattern to exclude (regex string), returns pandas dataframe with sample_id column and one (fa for fasta) or two (fq1 fq2, for fastq) file path columns. 
    """
    file_series = pd.Series(file_lst, dtype="str") #to fascilitate filtering
    ss_df = pd.DataFrame(dtype="str") #to store sample sheet
    assert mode in [0,1], f"Accepted mode values are 0 for fasta and 1 for fastq: {mode} was given."

    if mode == 1:  #If function is used to produce sample sheet from fasta files
        id_extractor = lambda x: os.path.basename(x).replace(generic_str, "") #extract id from string by replacing generic part
        id_series = file_series.apply(id_extractor) 
        if regex_str is not None: 
            id_series = id_series[id_series.str.contains(regex_str)] #additional sample id filtering based on regex was requested
            assert len(id_series) > 0, 'After filtering sample ids using regex no sample ids left'
        path_series = file_series[file_series.str.contains("|".join(id_series))].reset_index(drop=True) #getting corresponding paths to fastq files
        ss_df['sample_id'], ss_df['fa'] = id_series, path_series #adding to sample sheet dataframe
        return ss_df
    
    id_extractor = lambda x: re.sub(generic_str,"",os.path.basename(x)) #extract id from string by using regex
    id_series = file_series.apply(id_extractor).drop_duplicates(keep = "first").sort_values().reset_index(drop=True)

    if regex_str is not None: #additional sample id filtering based on regex was requested
        id_series = id_series[id_series.str.contains(regex_str)]
        assert len(id_series) > 0, 'After filtering sample ids using regex no sample ids left'
    read_1_dict, read_2_dict = {}, {} #to use python mapping to ensure correspondance between id and path

    for id in id_series:
        read_files = file_series[file_series.str.contains(id)].reset_index(drop=True).sort_values().reset_index(drop=True) #extract read paths
        read_1_dict[id] = read_files[0]
        read_2_dict[id] = read_files[1]
    ss_df['sample_id'] = id_series #adding to sample sheet dataframe
    ss_df['fq1'] = ss_df['sample_id'].map(read_1_dict)
    ss_df['fq2'] = ss_df['sample_id'].map(read_2_dict)

    return ss_df


def edit_sample_sheet(ss_df, info_dict, col_name):
    """
    Given sample sheet as pandas dataframe (object), a dictionary where each sample id is matched with information to be added (dict, values to be added as one column),
    and a new column name (str), returns a pandas dataframe (object), that contains new column where new information is added to the corresponding sample id.
    """
    ss_df[col_name] = ss_df["sample_id"].map(info_dict)
    return ss_df


def check_module_output(file_list):
    """
    Given (list) of paths to expected module output files, returns a dictionary where each file path is matched with the boolean (dict)
    indicating if it is present in the file system.
    """
    return {file: os.path.isfile(file) for file in file_list}
            

def read_config(config_path):
    """
    Given path to a config.yaml file, return a dictionary (dict) form of the yaml file.
    """
    with open(os.path.abspath(config_path), 'r') as yaml_handle:
        config_dict=yaml.safe_load(yaml_handle)
    return config_dict


def edit_config(config_dict, param, new_value):
    """
    Given a dictionary (dict) that is generated from config yaml file, a parameter that needs to be changes 
    and a new value of the parameter (string), return edited dictionary were the value of specified parameter is changed.
    (Adjusted from here: https://localcoder.org/recursively-replace-dictionary-values-with-matching-key)
    Return 0 if key was found and value changed, 1 otherwise.
    """
    if param in config_dict:
        config_dict[param] = new_value
        return 0 #this return is reached if key was found and value was changed
    
    for param, value in config_dict.items():
        if isinstance(value, dict):
            edit_config(value, param, new_value)
    return 1 #this return is reached only when all recursive calls are made and key is not found


def get_all_keys(input_dict, key_set=set()):
    """
    Given a nested dictionary (dict), return a (set) of all keys in that dictionary. 
    When called multiple times without passing new set object, all keys get saved to the same set.
    Default function call example: 
        get_all_keys(input_dict, set()) - returns all keys in input_dict, overwrites the content of key_set, if function was called before.
    Accumulating function call example: 
        get_all_keys(input_dict):
            if function is called for the first time - default behavior
            if function called multiple times - adds all keys in input_dict to the content of key_set 
            (assign key_set value to a variable when calling a function in order not to lose the reference to it)
    """
    for key, value in input_dict.items(): #start at the top level
        key_set.add(key) #add top keys
        if isinstance(value, dict): #if nested subdict detected
            get_all_keys(value, key_set) #make a recursive call
    return key_set #return is reached only when there are no recursive calls, hence all nested structure was parsed

    
def validate_config(config_dict, template_config_path='./config_modular.yaml'):
    """
    Given a dictionary (dict), return 0 if the structure of the dictionary corresponds to the config template structure (read from file),
    return 1 if some keys are missing in the dictionary, return 2 if some new keys are found in the dictionary.
    """
    template_config_file = read_config(template_config_path)
    valid_key_dict = {key:0 for key in get_all_keys(template_config_file, set())}
    found_keys = list(get_all_keys(config_dict, set()))
    for key in found_keys:
        if key not in valid_key_dict:
            return 2
        else:
            valid_key_dict[key] += 1
    if all(valid_key_dict.values()):
        return 0
    else:
        return 1


def write_config(config_dict, config_path):
    """
    Given a dictionary (dict) and a path to the new config file (str) write the contents to the new config file.
    """
    with open(config_path, "w+") as config_handle:
        yaml.dump(config_dict,config_handle)


def submit_module_job(module_name, config_path, output_dir, jobscript_path='./ardetype_jobscript.sh'):
    """
    Given snakemake module name (str), path to the job_script and path to the config file 
    perform job submition to RTU HPC cluster, returning bytestring, representing job id.
    """
    modules = {
        "core":os.path.abspath("./snakefiles/bact_core"),
        "shell":os.path.abspath("./snakefiles/bact_shell"),
        "tip":os.path.abspath("./snakefiles/bact_tip"),
        "shape":os.path.abspath("./snakefiles/bact_shape")
    }
    shutil.copy(jobscript_path, f'{output_dir}ardetype_jobscript.sh')

    try:
        job_id = subprocess.check_output(['qsub', '-F', f'{modules[module_name]} {config_path}', f'{output_dir}ardetype_jobscript.sh'])
    except subprocess.CalledProcessError as msg:
        sys.exit(f"Job submission error: {msg}")
    os.system(f"rm {output_dir}ardetype_jobscript.sh")
    return job_id


def check_job_completion(job_id, module_name, job_name="ardetype", sleeping_time=150, output_dir=None):
    """
    Given job id (bytestring) and sleeping time (in seconds) between checks (int), module name (str) and job name (str),
    checks if the job is complete every n seconds. Waiting is finished when the job status is C (complete). Optionally, 
    path to output directory (str) may be given to move the file with the job standard output to it.
    """
    job_id = job_id.decode('UTF-8').strip()
    search_string = f"{job_id}.*{job_name}.*{getuser()}.*[RCQ]"
    print(f'Going to sleep until ardetype/{module_name}/{job_id} job is finished')
    while True:
        qstat = os.popen("qstat").read()
        check_job = re.search(search_string,qstat).group(0)
        print(f"{check_job} : {time.ctime(time.time())}")
        if check_job[-1] == "C":
            print(f'Finished waiting: ardetype/{module_name}/{job_id} is complete')
            break
        time.sleep(sleeping_time)
    if output_dir is not None:
        job_report = f"{job_name}.o{job_id.split('.')[0]}"
        os.system(f"mv {job_report} {output_dir}/{module_name}_{job_name}_{job_id}.txt")


def install_snakemake():
    '''Function is used as a wrapper for bash script that checks if snakemake is installed and installs if absent.'''
    os.system(
    '''
    eval "$(conda shell.bash hook)"
    DEFAULT_ENV=/mnt/home/$(whoami)/.conda/envs/mamba_env/envs/snakemake
    SEARCH_SNAKEMAKE=$(conda env list | grep ${DEFAULT_ENV})
    if [ ${SEARCH_SNAKEMAKE} -ef ${DEFAULT_ENV} ]; then
        echo Running with --install_snakemake flag: Snakemake is already installed for this user
    else
        echo Running with --install_snakemake flag:
        conda create -n mamba_env -c conda-forge mamba
        conda activate mamba_env
        mamba create -c conda-forge -c bioconda -n snakemake snakemake
        conda activate snakemake
    fi    
    '''
    )


def run_module_cluster(module_name, config_path, cluster_config, job_count):
    '''
    Given snakemake module name (str) and path to the config file (str), path to the cluster config file (str)
    and the number of jobs to run in parallel (int) runs module on the login node of the HPC cluster, 
    allowing the snakemake to do job submissions to the computing nodes automatically.     
    '''
    modules = {
        "core":os.path.abspath("./snakefiles/bact_core"),
        "shell":os.path.abspath("./snakefiles/bact_shell"),
        "tip":os.path.abspath("./snakefiles/bact_tip"),
        "shape":os.path.abspath("./snakefiles/bact_shape")
    }
    qsub_command = '"qsub -N {cluster.jobname} -l nodes={cluster.nodes}:ppn={cluster.ppn},pmem={cluster.pmem},walltime={cluster.walltime} -q {cluster.queue} -j {cluster.jobout} -o {cluster.outdir} -V"'
    shell_command = f'''
    eval "$(conda shell.bash hook)";
    conda activate /mnt/home/$(whoami)/.conda/envs/mamba_env/envs/snakemake; 
    snakemake --jobs {job_count} --cluster-config {cluster_config} --cluster-cancel qdel --configfile {config_path} --snakefile {modules[module_name]} --keep-going --use-envmodules --use-conda --conda-frontend conda --rerun-incomplete --latency-wait 30 --cluster {qsub_command} --forceall'''
    try:
        subprocess.check_call(shell_command, shell=True)
    except subprocess.CalledProcessError as msg:
        sys.exit(f"Module process running error: {msg}")


def remove_invalid_samples(sample_sheet, module_name, output_dir):
    '''
    Given sample_sheet (DataFrame) and module_name (str), overwrites sample sheet in the specified output_dir (str), 
    removing samples that lack files, required by the module. If all samples are removed, returns None (float).
    '''
    module_requests = {
        "shell":"_contigs.fasta:True"
    }

    sample_sheet = sample_sheet[sample_sheet['check_note_core'].str.contains(module_requests[module_name])]
    if sample_sheet.empty:
        return 1
    else:
        sample_sheet.to_csv(f"{output_dir}sample_sheet.csv", header=True, index=False)


def parse_arguments():
    """
    Parse pre-defined set of arguments from the command line, returning a namespace (object),
    that allows accessing arguments as instance variables of namespace by their full name.
    """    
    ###Parsers
    parser = argparse.ArgumentParser(description='This is a wrapper script of ARDETYPE pipeline.', formatter_class=argparse.RawTextHelpFormatter)
    req_arg_grp = parser.add_argument_group('required arguments') #to display argument under required header in help message
    
    ###generic arguments
    #####Required
    req_arg_grp.add_argument('-m', '--mode',
        metavar='\b',
        help = """Selecting mode that allows to run specific modules of the pipeline:
        all - run all modules (starting from fastq files) (not active)
        core - run only bact_core module (starting from fastq files) 
        shell - run only bact_shell module (starting from fasta file) (not active)
        shell_tip - run bact_shell and bact_tip modules (starting from fasta file) (not active)
         """,
        default=None,
        required=True)

    #####Optional
    parser.add_argument('-c', '--config', metavar='\b', help = 'Path to the config file (if not supplied, the copy of the template config_modular.yaml will be used)', default="./config_modular.yaml", required=False)
    parser.add_argument('-r', '--skip_reporting', help = 'Use this flag to skip reporting trough bact_shape module (which will run by-default with any other option)', action='store_true')
    parser.add_argument('-o', '--output_dir', metavar='\b', help = 'Path to the output directory where the results will be stored (ardetype_output/ by-default).', default="ardetype_output/", required=False)
    parser.add_argument('-s', '--install_snakemake', help = 'Use this flag to install mamba and snakemake for the current HPC user, if it is not already installed.', action='store_true')
    parser.add_argument('-j', '--module_jobs', help='Use this flag to run modules as individual jobs on HPC (without submitting subjobs to computational nodes)', action='store_true')
    ###bact_core arguments

    #####Required
    req_arg_grp.add_argument('-i', '--input', metavar='\b', help = 'Path to directory that contains files to be analysed (all files in subdirectories are included).', default=None, required=True)
    #####Optional

    ###If no command-line arguments provided - display help and stop script excecution
    if len(sys.argv)==1: 
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args