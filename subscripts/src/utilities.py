import os, sys, yaml, pandas as pd, re, argparse, json, base64, requests, numpy as np, urllib, pandas as pd, concurrent.futures
from dateutil.relativedelta import relativedelta
from Bio import SeqIO, Entrez
from datetime import datetime
from pathlib import Path
from shutil import move

class Housekeeper:
    '''Class to contain methods that perform general housekeeping tasks for the pipeline, 
    like reading from/writing to certain file types, generating sample sheets etc.'''

    @staticmethod
    def parse_folder(folder_pth_str:str, file_fmt_str:str, substr_lst:list=None, regstr_lst:list=None):
        '''
        Given path to the folder (folder_pth_str) and file format (file_fmt_str), returns a list, 
        containing absolute paths to all files of specified format found in folder and subfolders,
        except for files that contain patterns to exclude (specified in regstr_lst) or substrings to exclude (specified in substr_lst).    
        '''
        name_series = pd.Series(dtype="str") #initialize pandas series to store path values
        for (root,_,files) in os.walk(folder_pth_str, topdown=True): #get list of file paths (from parent dir & subdir)
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

    @staticmethod
    def create_sample_sheet(file_lst:list, generic_str:str, regex_str:str=None, mode:int=0):
        """
        Given (list) of paths to files and a generic part of the file name (e.g. _contigs.fasta or _R[1,2]_001.fastq.gz string, regex expected for fastq), mode value (int 1 for fasta, 0 (default) for fastq)
        and a sample_id regex pattern to exclude (regex string), returns pandas dataframe with sample_id column and one (fa for fasta) or two (fq1 fq2, for fastq) file path columns. 
        """
        file_series = pd.Series(file_lst, dtype="str") #to facilitate filtering
        ss_df = pd.DataFrame(dtype="str") #to store sample sheet
        if mode not in [0,1]:
            raise Exception(f"utilities/create_sample_sheet: Accepted mode values are 0 for fasta and 1 for fastq: {mode} was given.") 

        if mode == 1:  #If function is used to produce sample sheet from fasta files
            id_extractor = lambda x: os.path.basename(x).replace(generic_str, "") #extract id from string by replacing generic part
            id_series = file_series.apply(id_extractor) 
            if regex_str is not None: 
                id_series = id_series[id_series.str.contains(regex_str)] #additional sample id filtering based on regex was requested
                if len(id_series) == 0:
                    raise Exception('utilities/create_sample_sheet: After filtering sample ids using regex no sample ids left')
            path_series = file_series[file_series.str.contains("|".join(id_series))].reset_index(drop=True) #getting corresponding paths to fastq files
            ss_df['sample_id'], ss_df['fa'] = id_series, path_series #adding to sample sheet dataframe
            return ss_df
        
        id_extractor = lambda x: re.sub(generic_str,"",os.path.basename(x)) #extract id from string by using regex
        id_series = file_series.apply(id_extractor).drop_duplicates(keep = "first").sort_values().reset_index(drop=True)
        if regex_str is not None: #additional sample id filtering based on regex was requested
            id_series = id_series[id_series.str.contains(regex_str)]
            if len(id_series) == 0:
                raise Exception('utilities/create_sample_sheet: After filtering sample ids using regex no sample ids left')
        read_1_dict, read_2_dict = {}, {} #to use python mapping to ensure correspondance between id and path

        for id in id_series:
            read_files = file_series[file_series.str.contains(id)].reset_index(drop=True).sort_values().reset_index(drop=True) #extract read paths
            read_1_dict[id] = read_files[0]
            read_2_dict[id] = read_files[1]
            
        ss_df['sample_id'] = id_series #adding to sample sheet dataframe
        ss_df['fq1'] = ss_df['sample_id'].map(read_1_dict)
        ss_df['fq2'] = ss_df['sample_id'].map(read_2_dict)
        return ss_df

    @staticmethod
    def map_new_column(ss_df:pd.DataFrame, info_dict:dict, id_column:str, new_col_name:str):
        """
        Given a pandas dataframe (object), a dictionary where each row in id_column is matched with information to be added (dict, values to be added as one column),
        and a new column name (str), returns a pandas dataframe (object), that contains new column where new information is added to the corresponding row of id_column.
        """
        if not isinstance(ss_df, pd.DataFrame): raise TypeError('Expected pandas.DataFrame as ss_df')
        elif not isinstance(info_dict, dict): raise TypeError('Expected dictionary as info_dict')
        elif id_column not in ss_df.columns: raise KeyError('id_column should be present in ss_df')
        elif not set(info_dict.keys()).intersection(set(ss_df[id_column])): raise KeyError('No overlap between ids in ss_df.id_column and info_dict')

        ss_df[new_col_name] = ss_df[id_column].map(info_dict)
        return ss_df

    @staticmethod
    def check_file_existance(file_list:list):
        """
        Given (list) of paths files, returns a dictionary where each file path is matched with the boolean (dict)
        indicating if it is present in the file system.
        """
        return {file: os.path.isfile(file) for file in file_list}
                
    @staticmethod
    def read_yaml(yaml_path:str):
        """
        Given path to a yaml file, returns a dictionary (dict) form of the yaml file.
        """
        with open(os.path.abspath(yaml_path), 'r') as yaml_handle:
            yaml_dict=yaml.safe_load(yaml_handle)
        return yaml_dict

    @staticmethod
    def edit_nested_dict(config_dict:dict, param, new_value):
        """
        Given a nested dictionary (dict), a parameter (key) that needs to be changed,
        and a new value of the parameter, returns edited dictionary were the value of specified parameter is changed.
        (Adjusted from here: https://localcoder.org/recursively-replace-dictionary-values-with-matching-key)
        Return 0 if key was found and value changed, None otherwise.
        """
        if param in config_dict:
            config_dict[param] = new_value
            return 0 #this return is reached if key was found and value was changed
        else:
            for value in config_dict.values():
                if isinstance(value, dict):
                    return Housekeeper.edit_nested_dict(value, param, new_value)

    @staticmethod
    def find_in_nested_dict(nested_dict:dict, key_sequence:list):
        '''
        Given a dictionary and an ordered sequence of keys in a form of list, returns value mapped to last key in sequence, by parsing the dictionary. 
        Raises exceptions if key is not found or non-dict value reached before last key in sequence is reached.
        '''
        if not isinstance(nested_dict,dict):
            raise TypeError('nested_dict should be a python dictionary')
        elif not hasattr(key_sequence, 'pop'):
            raise TypeError('key_sequence should have pop method defined')

        key = key_sequence.pop(0)
        try:
            if isinstance(nested_dict[key], dict) and len(key_sequence) != 0:
                tmp_dict = nested_dict[key]
            elif len(key_sequence) == 0:
                return nested_dict[key]
            elif not isinstance(nested_dict[key], dict):
                raise Exception('Problem with keys: reached non-dict value before processing all keys in sequence.')
        except KeyError:
            raise Exception(f'Problem with keys: {key} not found in nested_dict.')

        for key in key_sequence:
            try:
                if isinstance(tmp_dict[key], dict) and key != key_sequence[-1]:
                    tmp_dict = tmp_dict[key]
                elif key == key_sequence[-1]:
                    return tmp_dict[key]
                elif not isinstance(tmp_dict[key], dict):
                    raise LookupError('Problem with keys: reached non-dict value before processing all keys in sequence.')
            except KeyError:
                raise LookupError(f'Problem with keys: {key} not found in nested_dict.')

            
    @staticmethod
    def get_all_keys(input_dict:dict, key_set=set()):
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
                Housekeeper.get_all_keys(value, key_set) #make a recursive call
        return key_set #return is reached only when there are no recursive calls, hence all nested structure was parsed

    @staticmethod
    def validate_yaml(input_dict:dict, template_yaml_path:str='./config_files/yaml/config_modular.yaml'):
        """
        Given a dictionary (dict), return 0 if the structure of the dictionary corresponds to the yaml template structure (read from file),
        return 1 if some keys are missing in the dictionary, return 2 if some new keys are found in the dictionary.
        """
        template_yaml_path = Housekeeper.read_yaml(template_yaml_path)
        valid_key_dict = {key:0 for key in Housekeeper.get_all_keys(template_yaml_path, set())} #initializing to add 1 to all valid keys
        found_keys = list(Housekeeper.get_all_keys(input_dict, set())) #get list of all keys
        for key in found_keys:
            if key not in valid_key_dict: #if undefined key is found
                return 2
            else:   
                valid_key_dict[key] += 1 #if valid key is found - set its check value to 1 (True)
        if all(valid_key_dict.values()): #if all defined keys are found
            return 0
        else:   #if some defined keys are missing
            return 1

    @staticmethod
    def write_yaml(input_dict:dict, yaml_path:str):
        """
        Given a dictionary (dict) and a path to the new config file (str) write the contents to the new config file.
        """
        with open(yaml_path, "w+") as yaml_handle:
            yaml.dump(input_dict,yaml_handle)

    @staticmethod
    def write_json(input_dict:dict, json_path:str):
        """
        Given a dictionary (dict) and a path to the json file (str) writes the contents to file.
        """
        with open(json_path, "w+") as json_handle:
            json.dump(input_dict,json_handle)

    @staticmethod
    def install_snakemake():
        '''Function is used as a wrapper for bash script that checks if snakemake is installed and installs if absent.'''
        os.system(
        '''
        eval "$(conda shell.bash hook)"
        DEFAULT_ENV=/mnt/home/$(whoami)/.conda/envs/mamba_env/envs/snakemake$
        SEARCH_SNAKEMAKE=$(conda env list | grep -oP "${DEFAULT_ENV}")
        if [ ${SEARCH_SNAKEMAKE} -ef ${DEFAULT_ENV::-1} ]; then
            echo Running with --install_snakemake flag: Snakemake is already installed for this user
        else
            echo Running with --install_snakemake flag:
            conda create -n mamba_env
            conda activate mamba_env
            conda install python=3.9
            conda install -c conda-forge mamba
            mamba create -c conda-forge -c bioconda -n snakemake snakemake=7.6.1
            conda activate /mnt/home/$(whoami)/.conda/envs/mamba_env/envs/snakemake
            pip install PyYAML bs4 lxml
        fi    
        '''
        )

    @staticmethod
    def read_json_dict(json_path:str):
        '''
        Given path to a json file, returns python (dict) object.
        '''
        with open(json_path) as json_file:
            return json.load(json_file)

    @staticmethod
    def type_contigs_api(contigs_path:str, organism:str, scheme_num:int=0):
        '''
        Given path to a fasta file and the full name of the organism, sends POST request to the corresponding API.
        If returned status code is valid, returns dictionary, containing response details, else returns dictionary with status code and corresponding text.
        If no api is available for given organism returns 2. Optionally, scheme_num value might be given if more than one typing scheme is in use for the same organism.
        '''
        #this dictionary is the main interface for typing option selection - one organism - one option
        organisms = {
            'Listeria monocytogenes': [
                "https://bigsdb.pasteur.fr/api/db/pubmlst_listeria_seqdef/schemes/2/sequence",
                "https://bigsdb.pasteur.fr/api/db/pubmlst_listeria_seqdef/schemes/3/sequence"
                ],
            'Neisseria gonorrhoeae': [
                "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/71/sequence",
                "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/62/sequence"
                ],
            'Acinetobacter baumanii': [
                "https://rest.pubmlst.org/db/pubmlst_abaumannii_seqdef/schemes/1/sequence",
                "https://rest.pubmlst.org/db/pubmlst_abaumannii_seqdef/schemes/3/sequence"
                ],
            "Neisseria meningitidis":
                ["https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/47/sequence"],
            "Klebsiella pneumoniae":
                ["https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/schemes/15/sequence"],
            "Staphylococcus aureus":
                ["https://rest.pubmlst.org/db/pubmlst_saureus_seqdef/schemes/2/sequence"],
            "Salmonella enterica":
                ["https://rest.pubmlst.org/db/pubmlst_salmonella_seqdef/schemes/4/sequence"],
            "Escherichia coli":
                ["https://rest.pubmlst.org/db/pubmlst_escherichia_seqdef/schemes/6/sequence"],
            "Streptococcus pneumoniae":
                ["https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/schemes/2/sequence"]
        }
        try:
            url = organisms[organism][scheme_num] #trying to get url for the organism requested by the user
        except KeyError:
            return 2 #if url is not defined
        with open(contigs_path, 'r') as x: #if url is obtained
            fasta = x.read()
        req_start = '{"base64":true,"details":true,"sequence":"' #request wrapper as defined in API example
        req_stop = '"}'
        payload = req_start + base64.b64encode(fasta.encode()).decode() + req_stop #request + contig sequences
        response = requests.post(url, data=payload) #sending request to the API and saving responce to a variable
        if response.status_code == requests.codes.ok:   #if responce status code is valid (e.g. not 404)
            return response.json()  #return dictionary with the response contents
        else: #if return code indicates failure
            return {"status_code": response.status_code, "text": response.text} #return dictionary with corresponding status code

    @staticmethod
    def filter_contigs_length(input_multifasta_path:str, output_multifasta_path:str, minlen:int=500):
        '''
        Given path to multifasta, filters out contigs that are less than specified length (default 500 bp)
        and saves filtered contigs to a new multifasta file.
        '''
        SeqIO.write([record for record in SeqIO.parse(input_multifasta_path, "fasta") if len(record.seq) > minlen], output_multifasta_path, "fasta")

    @staticmethod
    def check_file_multiplicity(file_path_list:list):
        '''
        Function checks if list of file paths contains files that are multiplicated files (e.g. paired fastq files).
        Returns integer, indicating multiplicity of files in the list. Raises an error if multiplicity is greater than 1000.
        Assumes the same multiplicity for all files in the list (e.g. list sould only contain path to pair-end or single-end fastq, but not both)
        '''
        get_file_names = np.vectorize(lambda x: os.path.basename(x))
        file_names = get_file_names(np.array(file_path_list))
        current_array = np.array(list(file_names[0]))
        current_multiplicity, current_order = 1, 0
        for name in file_names[1:]:
            array = np.array(list(name))
            try:
                comparison_arr = array == current_array
                if comparison_arr.sum() == len(current_array):
                    return current_multiplicity

                comparison = comparison_arr.sum() >= len(current_array) - 3
                if comparison:
                    diff_index = np.where(~comparison_arr)[0][0]
                    diff_num = abs(int(current_array[diff_index]) - int(array[diff_index])) == 1
                    if diff_num:
                        current_multiplicity += 1
                        current_array = array                   
                    else:
                        return current_multiplicity
                else:
                    return current_multiplicity
            except ValueError:
                if int("".join(array[diff_index:diff_index+current_order+1])) // 10 == 0:
                    current_order += 1
                    current_multiplicity += 1
                    current_array = array
                else:
                    return current_multiplicity
        return current_multiplicity
        
    @staticmethod
    def parse_arguments(arg_dict:dict):
        """
        Parse pre-defined set of arguments from the command line, returning a namespace (object),
        that allows accessing arguments as instance variables of namespace by their full name. 
        Expected arg_dict structure:
        {
            'description':'',
            'required_arguments':[
                ['--r1','--first_required','r1_help message'],
                ...
            ],
            'optional':{
                'arguments':[
                    ['--o1','--first_optional','o1_help message', 'o1_default_value'],
                    ...
                ],
                'flags':[
                    ['--f1','--first_flag','f1_help message'],
                    ...
                ]
            }
        }
        """

        #Argument parsers
        parser = argparse.ArgumentParser(description=arg_dict['description'], formatter_class=argparse.RawTextHelpFormatter)
        req_arg_grp = parser.add_argument_group('Required arguments') #to display argument under required header in help message
        
        #Required
        for req_arg in arg_dict['required_arguments']: 
            req_arg_grp.add_argument(req_arg[0], req_arg[1], metavar='\b', help=req_arg[2],default=None, required=True)

        #Optional
        #Arguments
        for opt_arg in arg_dict['optional']['arguments']: 
            parser.add_argument(opt_arg[0], opt_arg[1], metavar='\b', help=opt_arg[2], default=opt_arg[3], required=False)
        
        #Flags
        for flag in arg_dict['optional']['flags']: parser.add_argument(flag[0], flag[1], help = flag[2], action='store_true')

        #If no arguments provided - display help and stop the script
        if len(sys.argv)==1: 
            parser.print_help(sys.stderr)
            sys.exit(1)
        args = parser.parse_args()
        return args

    @staticmethod
    def asign_perm_rec(path_to_folder:str, linux_permissions:str="775"):
        '''
        Given path_to_folder string, representing a path in linux-based system,
        recursively assigns permissions (775 by-default) to all files in the folder.
        '''
        for root, dirs, files in os.walk(path_to_folder):
            for d in dirs:
                try:
                    os.chmod(os.path.join(root, d), int(linux_permissions, 8))
                except PermissionError:
                    continue
            for f in files:
                try:
                    os.chmod(os.path.join(root, f), int(linux_permissions, 8))
                except PermissionError:
                    continue
        os.chmod(path_to_folder, int(linux_permissions, 8))

    @staticmethod
    def extract_log_id(path_to_log:str, pattern_to_search:str="wildcards: sample_id_pattern=.*"):
        '''Given path to log file and pattern to search, returns the search result or False if search failed.'''
        with open(path_to_log, "r+") as file: contents = file.read()
        found = re.search(pattern_to_search, contents)
        if found is not None:
            found = found[0].replace("wildcards: sample_id_pattern=", "")
            return found
        else:
            return False

    @staticmethod
    def rename_file(path_to_file:str, pattern_to_add:str, target_folder_path:str=f"{os.path.dirname(Path(__file__).parents[0].absolute())}/ardetype_job_logs"):
        '''
        Given path to file and pattern to add in file name, 
        saves file in the same directory with pattern added to its name.
        '''
        new_path = f"{target_folder_path}/{pattern_to_add}_{os.path.basename(path_to_file)}"
        move(path_to_file,new_path)


    @staticmethod 
    def name_job_logs(pipeline_name:str):
        '''Given path to the pipeline_name_job_logs folder, adds sample id to the name of each log file (skips the file if sample id is not found in it).'''
        path_to_log_dir:str=f"{os.path.dirname(Path(__file__).parents[1].absolute())}/{pipeline_name}_job_logs"
        log_list = [f'{path_to_log_dir}/{f}' for f in os.listdir(path_to_log_dir) if f"_{pipeline_name}" not in f and "_default" not in f] #if pipeline name in log name is preceeded by _ it is already annotated
        print(f'\nAdding sample ids to log file names.\n')
        log_count = len(log_list)
        for i,log in enumerate(log_list):
            Housekeeper.printProgressBar(i+1, log_count, prefix = 'Progress:', suffix = 'Complete', length = 50)
            try:
                sample_id = Housekeeper.extract_log_id(log)
            except PermissionError:
                continue
            if sample_id and sample_id not in log: 
                Housekeeper.rename_file(log, sample_id, path_to_log_dir) #if no sample id is found, the extract_log_id returns False



    @staticmethod
    def remove_old_files(path_to_folder:str, valid_days:int=30):
        '''
        Given path to folder, removes all files that were created more than 30 days ago.
        By providing number of days beyond which the files are discarded (valid_days), 
        removes all the files in folder that are older than that number of days.
        '''
        os.chdir(path_to_folder)
        file_list = os.listdir(path_to_folder)
        today=datetime.now() #current time
        print(f'\nRemoving log files older than {valid_days} days.\n')
        file_count = len(file_list)
        for i,file in enumerate(file_list):
            file_date = datetime.fromtimestamp(os.path.getctime(f"./{file}")) #getting creation time of a file
            dif_days=(today-file_date).days #checking how many days passed since the file was created
            if dif_days > valid_days: #if file is older - delete it
                os.remove(f'./{file}')
            Housekeeper.printProgressBar(i+1, file_count, prefix = 'Progress:', suffix = 'Complete', length = 50)


    @staticmethod
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


    @staticmethod
    def find_job_logs(pipeline_name:str, logs_to_skip:list=[]) -> list:
        '''
        Given pipeline name, returns list of paths (as str) to all log files that contain pipeline name as substring in file name.
        Returns empty list if none is found. 
        '''
        path_to_log_dir=f"{os.path.dirname(Path(__file__).parents[1].absolute())}/{pipeline_name}_job_logs" #get path to log folder - static for default pipeline template
        processed_log_set = set(logs_to_skip) #to use set operations for speedup
        path_joiner = lambda p: os.path.join(path_to_log_dir, p) #helper function to apply map instead of using for loop
        full_log_set = set(map(path_joiner, os.listdir(path_to_log_dir))) #applying helper to all log paths to get set of full paths
        unprocessed_logs = list(full_log_set - processed_log_set) #using set operations to keep only paths to unprocessed logs
        return unprocessed_logs
                

    @staticmethod
    def parse_snakemake_log(path_to_log:str) -> pd.DataFrame:
        '''Given path to a snakemake log, returns pandas dataframe containing the following fields:
            log_path <str>
            job_name <str>
            input_size_gb <int>
            is_failed <int8>
            mem_gb_snakemake <int>
            cpu_snakemake <int>
            mem_gb_req <int>
            time_sec_req <int>
            ctime_theor <float,2>
            start_time <datetime>
            time_sec_total <int>
            ctime_real <float,2>
            Eff <float,2>
        If parsing fails, returns empty pandas dataframe.
        '''

        ###Static variables
        smk_job_hdr = 'Building DAG of jobs...'
        smk_dt_fmt = '%b  %d %H:%M:%S %Y' #expected in snakemake log
        strt_dt_fmt = '%Y:%m:%d-%H:%M:%S' #stored in dataframe
        failed_job_tags = [
            'Will exit after finishing currently running jobs (scheduler).', #snakemake stopped
            '=>> PBS: job killed:', #resource limit exceeded
            'Exiting because a job execution failed. Look above for error message', #job failed
            'Error in rule',
            ]
        smk_dt_slc = [5,-1] #slice on a string
        strt_tm_idx = 0 #in content list
        end_tm_idx = -3 #in content list
        job_nm_idx = 1 #in content list
        jb_nm_slc = [5,-1] #slice on a string
        cntt_strt_idx = 6 #in content list


        ###Reading log file into list
        with open(path_to_log, 'r+') as handle:
            header = handle.readline().strip() #Read 1st line without separator characters
            if header == smk_job_hdr: #If 1st line is valid
                contents = list(map(str.strip, handle.readlines()))[cntt_strt_idx:] #Reading the remaining contents and removing separator characters
            else:
                return pd.DataFrame() #Log is not valid - return empty df

        ###Parsing data from cluster config file
        cluster_config_dict = Housekeeper.read_yaml(f"{os.path.dirname(Path(__file__).parents[1].absolute())}/config_files/yaml/cluster.yaml")
        try:
            ###Pre-processing contents
            is_failed = 1 if any([tag in row for tag in failed_job_tags for row in contents]) else 0
            job_start_time = datetime.strptime(contents[strt_tm_idx][smk_dt_slc[0]:smk_dt_slc[1]], smk_dt_fmt) #Parsing start time into datetime object
            job_end_time = datetime.strptime(contents[end_tm_idx][smk_dt_slc[0]:smk_dt_slc[1]], smk_dt_fmt) if not is_failed else datetime.fromtimestamp(os.path.getctime(path_to_log)) #Parsing end time into datetime object
            time_sec_total = relativedelta(job_end_time, job_start_time)
            time_sec_total = time_sec_total.hours*3600+time_sec_total.minutes*60+time_sec_total.seconds #getting job real runtime in seconds according to snakemake
            job_start_time = datetime.strftime(job_start_time, strt_dt_fmt) #converting start time to string to store in a dataframe
            job_name = contents[job_nm_idx][jb_nm_slc[0]:jb_nm_slc[1]] #Getting job name
            smk_ram = int(contents[6].replace('resources: ',"").split(", ")[0].replace('mem_mb=',''))/1024 if 'resources:' in contents[6] else int(contents[7].replace('resources: ',"").split(",")[0].replace('mem_mb=',''))/1024
            smk_threads = int(contents[6].replace('threads: ',"")) if 'threads' in contents[6] else 1
            sample_id = contents[5].replace("wildcards: sample_id_pattern=","")

            ###Getting data from cluster config
            procs_req = int(cluster_config_dict[job_name]['procs'])
            mem_gb_req = int(cluster_config_dict[job_name]['pmem'].replace('mb',''))*procs_req/1024 if "pmem" in cluster_config_dict[job_name] else int(cluster_config_dict[job_name]['mem'].replace('mb',''))/1024
            time_sec_req = relativedelta(hours=int(
                cluster_config_dict[job_name]['walltime'].split(':')[0]),
                minutes=int(cluster_config_dict[job_name]['walltime'].split(':')[1]),
                seconds=int(cluster_config_dict[job_name]['walltime'].split(':')[2]))
            time_sec_req = time_sec_req.hours*3600+time_sec_req.minutes*60+time_sec_req.seconds
            Eff = round(100*time_sec_total/time_sec_req,2)
            if sample_id not in path_to_log: path_to_log = path_to_log.replace('_job_logs/', f'_job_logs/{sample_id}_')
            
            df = pd.DataFrame({
                "log_path":[path_to_log],
                "job_name":[job_name],
                "sample_id":[sample_id],
                "is_failed":[is_failed],
                "mem_gb_snakemake":[smk_ram],
                "cpu_snakemake":[smk_threads],
                "mem_gb_req":[mem_gb_req],
                "time_sec_req":[time_sec_req],
                "start_time":[job_start_time],
                "time_sec_total":[time_sec_total],
                "Eff":[Eff]
            })
            return df
        except:
            return pd.DataFrame()
        

    @staticmethod
    def aggregate_job_logs(log_path_list:list, procs:int=24) -> pd.DataFrame:
        '''
        Given list of paths to log files of individual jobs of the pipeline, aggregates the log data in the pandas dataframe.
        Returns empty dataframe if input list is empty. Prints warnings for files that could not be properly parsed by the extractor function.
        '''
        aggr_df = pd.DataFrame({
                    "log_path":[],
                    "job_name":[],
                    "sample_id":[],
                    "is_failed":[],
                    "mem_gb_snakemake":[],
                    "cpu_snakemake":[],
                    "mem_gb_req":[],
                    "time_sec_req":[],
                    "time_sec_total":[],
                    "Eff":[]
                })
        with concurrent.futures.ProcessPoolExecutor(max_workers=procs) as executor:
            results = [executor.submit(Housekeeper.parse_snakemake_log, file_path) for file_path in log_path_list]
            processed_count = 0 #TO VIEW PROGRESS
            print(f'Aggregating job log information: ')
            for output in concurrent.futures.as_completed(results):
                try:
                    if not output.result().empty:
                        aggr_df = pd.concat([aggr_df, output.result()], sort=False)
                except PermissionError:
                    pass
                processed_count += 1 #counting processed files
                Housekeeper.printProgressBar(processed_count, len(log_path_list), prefix = 'Progress:', suffix = 'Complete', length = 50)

        return aggr_df


    @staticmethod
    def update_log_summary(notebook_path:str, env_path:str, output_dir:str) -> None:
        '''Running the notebook from the command line with nbconvert.'''
        os.system(
            f'''
            eval "$(conda shell.bash hook)" 
            conda activate {env_path}
            jupyter nbconvert --to html --execute --no-input {notebook_path} --output-dir={output_dir}
            chmod 775 {output_dir} 2>/dev/null
            ''')


    @staticmethod
    def update_log_history(pipeline_name:str) -> None:
        '''Updating current job log table & plots'''
        job_log_dir = f"{os.path.dirname(Path(__file__).parents[1].absolute())}/{pipeline_name}_job_logs"
        file_search = [path for path in os.listdir(f"{job_log_dir}") if '-log_aggregate_' in path]
        if file_search: current_file = file_search[0]
        else: current_file = 'none'
        if os.path.isfile(f'{job_log_dir}/{current_file}'): 
            current_df = pd.read_csv(f'{job_log_dir}/{current_file}')
        else:
            current_df = None
        path_list = Housekeeper.find_job_logs(pipeline_name, logs_to_skip=list(current_df['log_path']) if current_df is not None else [])
        new_log_df = Housekeeper.aggregate_job_logs(log_path_list=path_list)
        if current_df is not None and not new_log_df.empty: 
            updated_df = pd.concat([current_df, new_log_df], sort=False)
        elif new_log_df.empty:
            updated_df = current_df
        else:
            updated_df = new_log_df
        tstemp = datetime.now().strftime("%Y-%m-%d")
        new_file_name = f'{tstemp}-log_aggregate_{pipeline_name}.csv'
        updated_df.drop_duplicates(subset=['log_path'], keep='first', inplace=True)
        if file_search: os.remove(f'{job_log_dir}/{current_file}')
        updated_df.to_csv(f'{job_log_dir}/{new_file_name}', header=True, index=False)





class Query_ncbi:
    '''Class is set to contain methods to send requests to ncbi nucleotied database via biopython Entrez API.'''
    email = "jevgenijs.bodrenko@aslimnica.lv"
    database = 'nucleotide'

    @classmethod #uses email and database defined in Query
    def check_output(valid_accession:str):
        '''Method checks if record exists in NCBI database given accession number.'''
        Entrez.email = Query_ncbi.email
        try:
            handle = Entrez.efetch(db=Query_ncbi.database, id = valid_accession, rettype="acc") #Attempts to open connection to NCBI database
            data = handle.read() #Reads up-to-date accession from NCBI
            handle.close() #Closes connection
            if valid_accession in data: #Checks if found accession contains query
                return True
            else: 
                return False
        except urllib.error.HTTPError as e: #If query is invalid, the exception is thrown
            if str(e) == "HTTP Error 400: Bad Request":
                return False

    @classmethod
    def get_fasta(cls, valid_accession):
        '''Method attempts to get fasta sequence (header + sequence separately) given accession number.'''
        Entrez.email = Query_ncbi.email
        handle = Entrez.efetch(db=Query_ncbi.database, id = valid_accession, rettype="fasta") #Attempts to open connection to NCBI database
        record = SeqIO.read(handle, "fasta") #Reads sequence in fasta format from ncbi
        handle.close() #Closes connection
        return record.seq, record.id

    @classmethod
    def get_taxonomy(cls, valid_accession):
        '''Method attempts to get taxonomy information(list) from NCBI database given accession number.'''
        Entrez.email = cls.email
        handle = Entrez.efetch(db=Query_ncbi.database, id = valid_accession, rettype="gb") #Attempts to open connection to NCBI database
        tax_data = SeqIO.read(handle, format='genbank').annotations['taxonomy']
        handle.close() #Closes connection
        return tax_data