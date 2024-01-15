from .utilities import Housekeeper as hk
import os, warnings, re, subprocess, shutil, time, pandas as pd, glob, sys
from datetime import datetime
from itertools import chain
from getpass import getuser
from pathlib import Path
from shutil import move

#Suppressing pandas warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)


#Reading data used to build module objects
module_data   = hk.read_json_dict(f'{os.path.dirname(Path(__file__).parents[1].absolute())}/config_files/json/module_data.json')
pipeline_path = os.path.dirname(Path(__file__).parents[1].absolute())


####################
# Class definition
####################

class Module:
    '''Class represents single module of the ardetype pipeline'''

    def __init__(
            self, 
            module_name         : str, 
            input_path          : str, 
            module_config, 
            output_path         : str, 
            run_mode            : bool, 
            job_name            : str, 
            patterns            : dict, 
            targets             : list, 
            requests            : dict, 
            snakefile_path      : str, 
            cluster_config_path : str,
            retry_times         : int,
            dry_run             : bool, 
            force_all           : bool, 
            rule_graph          : bool, 
            pack_output         : bool, 
            rules_to_rerun      : list,
            unpack_output:bool) -> None:
        
        self.run_mode            = run_mode #If true, snakemake will be run as single job, else - will run as job submitter on the login node
        self.job_id              = None  #Will be added if self.run_mode is True and job was submitted to HPC; filled by submit_module_job
        self.taxonomy_dict       = None   #Required if module creates different targets for different samples based on taxonomy information; filled by add_taxonomy_column
        self.module_name         = module_name #To be used in configuration file & sample_sheet file + to connect between modules (using remove_invalid_samples)
        self.input_path          = input_path #to the folder containing fasta/fastq.gz files
        self.output_path         = f"{os.path.abspath(output_path)}/" #Path to the output folder, where files will be saved (converted to full path)
        self.target_list         = None #List of all target files the module expects to create; filled by fill_target_list
        self.sample_sheet        = None #to store current state of sample_sheet dataframe; filled by create_sample_sheet; altered by fill_sample_sheet & receive_sample_sheet
        self.aggr_taxonomy_path  = f'{os.path.abspath(self.output_path)}/{self.module_name}_aggregated_taxonomy.json' #where to look for top kraken2 hits if snakemake will produce it; used by add_taxonomy_column
        self.config_file_path    = f'{os.path.abspath(self.output_path)}/config.yaml' #where to look for operational copy of the configuration file; used by submit_module_job & run_module_cluster
        self.cluster_config_path = cluster_config_path #where to look for job resource definition file; used by run_module_cluster
        self.config_file         = hk.read_yaml(module_config) if isinstance(module_config, str) else module_config #read module configuration from file if string is supplied (path expected); else - reads dictionary; used by add_module_targets, add_output_dir, write_module_config
        self.input_dict          = {} #to store input file paths for each file extension; used by fill_input_dict, fill_sample_sheet, add_fasta_samples
        self.patterns            = patterns #to store file extension patterns of expected input files; used by fill_input_dict; fill_sample_sheet
        self.job_name            = job_name #to store job name if self.run_mode is True; used by check_job_completion
        self.targets             = targets #to store file extensions of expected output files; used by fill_target_list, check_module_output
        self.requests            = requests #to store file extensions for files that are neccessary to run the modules; used by remove_invalid_samples
        self.snakefile_path      = snakefile_path #to the rule file to be run as single job on HPC if self.run_mode is True; used by submit_module_job
        self.dry_run             = "-np" if dry_run else "" #to store dry-run flag if it is supplied, else empty string is stored
        self.retry_times         = retry_times #number of times snakemake will attempt to rerun failed jobs (default=3); used by run_module_cluster
        self.force_all           = "--forceall" if force_all else "" #to store forceall flag if it is supplied, else empty string is stored
        self.rule_graph          = f"--rulegraph | dot -Tpdf > {self.module_name}.pdf" if rule_graph else "" #to store rule_graph flag if it is supplied, else empty string is stored
        self.unpack_output       = True if force_all else unpack_output #used to move files outside sample folders and do a rerun; used by unfold_output
        self.removed_samples     = pd.DataFrame() #to store dataframe containing information about samples that were deemed invalid by the module
        self.pack_output         = pack_output #switch to control putting output files into one folder named after sample_id; used by fold_output
        self.cleanup_dict        = {} #to map origin paths of input files to path in working directory; filled by move_to_wd; used by clear_working_directory
        self.status_script       = f"{os.path.dirname(Path(__file__).parents[0].absolute())}/pbs-status.py"
        self.failed_stamp        = None #added if module has failed to produce requested files for 1 or more steps of the workflow
        self.rules_to_rerun      = [rule for rule in rules_to_rerun if rule in self._get_rule_names_from_snakefile(self.snakefile_path)] if rules_to_rerun is not None else []
        self.force_specific      = f"-R {' '.join(self.rules_to_rerun)}" if self.rules_to_rerun else ""

    @staticmethod
    def _get_rule_names_from_snakefile(snakefile_path):
        rule_names = []
        with open(snakefile_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('rule '):
                    rule_name = line.split(' ')[1].strip(':')
                    rule_names.append(rule_name)
        return rule_names


    def fill_input_dict(self, substring_list=['reads_unclassified', 'reads_classified'], mixed:bool=False, empty:bool=False):
        '''Fills self.input_dict using self.input_path and self.module_name by
        mapping each file format to the list of files of that format, found in the self.input_path, 
        excluding files that contain substrings in their names (supply None to avoid excluding files).
        If some files of required format are missing, raises an exception, indicating missing file format.'''
        if not mixed:
            for format in self.patterns['inputs']: 
                self.input_dict[format] = hk.parse_folder(self.input_path,substr_lst=substring_list, file_fmt_str=format)
                if not self.input_dict[format]: raise Exception(f'Missing {format} files in input directory')
        elif mixed or empty:
            for format in self.patterns['inputs']['required']: 
                self.input_dict[format] = hk.parse_folder(self.input_path,substr_lst=substring_list, file_fmt_str=format)
                if not self.input_dict[format]: raise Exception(f'Missing {format} files in input directory')
            if not empty:
                for format in self.patterns['inputs']['optional']:
                    parsed_files = hk.parse_folder(self.input_path,substr_lst=substring_list, file_fmt_str=format)
                    if parsed_files:
                        self.input_dict[format] = parsed_files

    def fill_sample_sheet(self):
        '''
        Initializes self.sample_sheet to pandas dataframe, using self.input_dict and self.module_name (restricted to fastq & fasta inputs).
        '''
        ###Development - create sample sheet from paired or unpaired files; update sample sheet with paired or unpaired files
        if len(self.input_dict) < 2: #only one file extension is used - assumed fastq.gz
            self.sample_sheet = hk.create_sample_sheet(self.input_dict["001.fastq.gz"],self.patterns['sample_sheet'],mode=0)
        else: #fastq & fasta assumed
            self.sample_sheet = hk.create_sample_sheet(self.input_dict["001.fastq.gz"],self.patterns['sample_sheet'],mode=0)
            fasta_dict = {re.sub("_contigs.fasta","",os.path.basename(contig)):contig for contig in self.input_dict["_contigs.fasta"]}
            self.sample_sheet = hk.map_new_column(self.sample_sheet,fasta_dict,'sample_id','fa')


    def fill_target_list(self, taxonomy_based:bool=False, mixed:bool=False, empty:bool=False):
        '''Fills self.target_list using data stored in self.sample_sheet instance variable.'''
        if taxonomy_based:#specific targets for each species
            self.target_list = [f'{self.output_path}{id}{tmpl}' for idx, id in enumerate(self.sample_sheet['sample_id']) for tmpl in self.targets[self.sample_sheet['taxonomy'][idx]]]
        elif mixed or empty:#both species-specific and non-specific targets
            self.target_list = [f'{self.output_path}{id}{tmpl}' for id in self.sample_sheet['sample_id'].to_list()+self.removed_samples['sample_id'].to_list() for tmpl in self.targets['general']]
            if not empty:
                self.target_list += [f'{self.output_path}{id}{tmpl}' for idx, id in enumerate(self.sample_sheet['sample_id']) for tmpl in self.targets[self.sample_sheet['taxonomy'][idx]]]
        else:#only non-specific targets
            self.target_list = [f'{self.output_path}{id}{tmpl}' for id in self.sample_sheet['sample_id'] for tmpl in self.targets]
            

    def make_output_dir(self):
        '''Creates output directory (if not present in the file system) using self.output_path.'''
        os.makedirs(self.output_path, exist_ok=True)


    def write_sample_sheet(self):
        '''Creates sample_sheet.csv file in the self.output_path folder, using self.sample_sheet.'''
        self.sample_sheet.to_csv(f"{self.output_path}sample_sheet.csv", header=True, index=False)


    def add_module_targets(self):
        '''Updates self.config_file, using self.module_name.'''
        output_code = hk.edit_nested_dict(config_dict=self.config_file, param=f"{self.module_name}_target_files", new_value=self.target_list)
        validation_code = hk.validate_yaml(self.config_file)
        if not output_code == 0: raise Exception(f'Config editing failed with error code {output_code}')
        elif not validation_code == 0: raise Exception(f'Config validation failed with error code {validation_code}')


    def add_output_dir(self):
        '''Updates self.config_file using self.output_path.'''
        output_code = hk.edit_nested_dict(config_dict=self.config_file, param="output_directory", new_value=self.output_path)
        validation_code = hk.validate_yaml(self.config_file)
        if not output_code == 0: raise Exception(f'Config editing failed with error code {output_code}')
        elif not validation_code == 0: raise Exception(f'Config validation failed with error code {validation_code}')


    def write_module_config(self):
        '''Writes self.config_file to the self.output_path'''
        hk.write_yaml(self.config_file, f'{self.output_path}config.yaml')


    def check_module_output(self, mixed:bool=False):
        '''Checks if output files are generated according to self.module_name and adds check_note_{self.module_name} column 
        to the self.sample_sheet dataframe, where boolean value is stored for each expected file.'''
        ###Development - Automatically scale dirs_up depending on input structure - currently two dirs up max
        check_dict = hk.check_file_existance(file_list=self.target_list)
        if mixed:
            id_check_dict = {id:"" for id in self.sample_sheet['sample_id'].to_list()+self.removed_samples['sample_id'].to_list()}
        else:
            id_check_dict = {id:"" for id in self.sample_sheet['sample_id']}
        for file in check_dict:
            two_dirs_up = os.path.basename(os.path.dirname(os.path.dirname(file)))+"/"+os.path.basename(os.path.dirname(file))+"/"+os.path.basename(file) #required for outputs where directory patterns are defined in addition to file extensions
            if isinstance(self.targets, list): #if only un-specific targets are supplied
                id = os.path.basename(re.sub("("+"|".join(self.targets)+")","",two_dirs_up))
            elif isinstance(self.targets, dict): #if taxonomy-based targets are supplied - all species-specific target lists are to be merged into one list using chain.from_iterables
                id = os.path.basename(re.sub("("+"|".join(chain.from_iterable(self.targets.values()))+")","",two_dirs_up))
            id_check_dict[id] += f"|{file}:{check_dict[file]}"
        self.sample_sheet = hk.map_new_column(self.sample_sheet, id_check_dict, 'sample_id', f"check_note_{self.module_name}")


    def supply_sample_sheet(self, removed:bool=False): #getter, may not be required now as all variables are public, but makes it easier to encapsulate later, if needed
        '''Returns self.sample_sheet dataframe object.'''
        return self.sample_sheet


    def receive_sample_sheet(self, sample_sheet:pd.DataFrame): #setter, may not be required now as all variables are public, but makes it easier to encapsulate later, if needed
        '''Inializes self.sample_sheet with external sample_sheet dataframe (used to connect modules).'''
        self.sample_sheet = sample_sheet


    def add_fasta_samples(self):
        '''Adds fa column with _contigs.fasta files to the self.sample_sheet dataframe.'''
        fasta_dict = {re.sub("_contigs.fasta","",os.path.basename(contig)):contig for contig in self.input_dict["_contigs.fasta"]}
        self.sample_sheet = hk.map_new_column(self.sample_sheet,fasta_dict,'sample_id','fa')


    def remove_invalid_samples(self, connect_from_module_name:str, taxonomy_only:bool=False):
        '''
        Removes samples that lack files, required by the current module, given supplier module name.
        If all samples are removed, returns 1 (int).
        '''
        if not taxonomy_only:
            if self.requests['check'] is not None: #if module requires certain file types to run rules that are not taxonomy-specific
                if isinstance(self.requests['check'],str): #if only one requirement is provided as string
                    self.removed_samples = self.sample_sheet[~self.sample_sheet[f'check_note_{connect_from_module_name}'].str.contains(self.requests['check'])].reset_index(drop=True)
                    self.sample_sheet = self.sample_sheet[self.sample_sheet[f'check_note_{connect_from_module_name}'].str.contains(self.requests['check'])].reset_index(drop=True)
                elif isinstance(self.requests['check'], list):
                    for request in self.requests['check']: #if many requirements are provided as list of stings
                        self.removed_samples = self.sample_sheet[~self.sample_sheet[f'check_note_{connect_from_module_name}'].str.contains(request)].reset_index(drop=True)
                        self.sample_sheet = self.sample_sheet[self.sample_sheet[f'check_note_{connect_from_module_name}'].str.contains(request)].reset_index(drop=True)
                elif isinstance(self.requests['check'], dict):
                    for request in self.requests['check']: #if many requirements are provided as list of stings
                        for check in self.requests['check'][connect_from_module_name]:
                            self.removed_samples = self.sample_sheet[~self.sample_sheet[f'check_note_{connect_from_module_name}'].str.contains(check)].reset_index(drop=True)
                            self.sample_sheet = self.sample_sheet[self.sample_sheet[f'check_note_{connect_from_module_name}'].str.contains(check)].reset_index(drop=True)
            if self.requests['taxonomy'] is not None:   #if module can run only certain taxonomy-specific rules (but not others) and accepted species are supplied as list of strings
                self.removed_samples = self.removed_samples.append(self.sample_sheet[~self.sample_sheet['taxonomy'].str.contains("("+"|".join(self.requests['taxonomy'])+")")].reset_index(drop=True))
                self.sample_sheet = self.sample_sheet[self.sample_sheet['taxonomy'].str.contains("("+"|".join(self.requests['taxonomy'])+")")].reset_index(drop=True)
            if self.sample_sheet.empty: #if all samples were removed by filtering - e.g. files are missing and/or no taxonomy-based rules are available for detected organisms
                return 1
        else:
            if self.requests['taxonomy'] is not None:   #if module can run only certain taxonomy-specific rules (but not others) and accepted species are supplied as list of strings
                self.removed_samples = self.removed_samples.append(self.sample_sheet[~self.sample_sheet['taxonomy'].str.contains("("+"|".join(self.requests['taxonomy'])+")")].reset_index(drop=True))
                self.sample_sheet = self.sample_sheet[self.sample_sheet['taxonomy'].str.contains("("+"|".join(self.requests['taxonomy'])+")")].reset_index(drop=True)
            if self.sample_sheet.empty: #if all samples were removed by filtering - e.g. files are missing and/or no taxonomy-based rules are available for detected organisms
                return 1


    def save_removed(self):
        '''Generates a csv file in self.output_path folder, containing information about samples that were
        filtered as invalid by the module (see remove_invalid_samples). Does nothing if self.removed_samples is empty.'''
        if not self.removed_samples.empty: 
            self.removed_samples.to_csv(f"{self.output_path}removed_samples_{self.module_name}.csv", header=True, index=False)
            return self.removed_samples
    

    def submit_module_job(self, jobscript_path):
        """
        Submits module as a job to HPC cluster, given path to the job_script, setting self.job_id to the bytestring representing job_id.
        All rules are run sequentially using same set of resources.
        """
        shutil.copy(jobscript_path, f'{self.output_path}ardetype_jobscript.sh') #to avoid running the template file

        try:
            self.job_id = subprocess.check_output(['qsub', '-F', f'{self.snakefile_path} {self.config_file_path}', f'{self.output_path}ardetype_jobscript.sh'])
            os.remove(f'{self.output_path}ardetype_jobscript.sh') #cleanup
            # os.system(f"rm {self.output_path}ardetype_jobscript.sh") #cleanup
        except subprocess.CalledProcessError as msg:
            raise Exception(f"{self.module_name} job submission error: {msg}")
        

    def check_job_completion(self, sleeping_time=5):
        """
        Checks if the job is complete every n seconds, given sleeping time (in seconds) between checks (int). 
        Waiting is finished when the job status is C (complete), then the file with job std(out/err) is moved to self.output_path.
        """
        self.job_id = self.job_id.decode('UTF-8').strip() #job id is returned as byte string
        search_string = f"{self.job_id}.*{getuser()}.*[RCQ]" #regex to check the current job status (Running, Complete or Queued)
        print(f'Going to sleep until ardetype/{self.module_name}/{self.job_id} job is finished') #Informing the user
        while True: #event-driven check
            qstat = os.popen("qstat").read() #read qstat output
            check_job = re.search(search_string,qstat).group(0) #check job status
            print(f"{check_job} : {time.ctime(time.time())}") #inform the user about last check time
            if check_job[-1] == "C": #if job is complete - stop waiting
                print(f'Finished waiting: ardetype/{self.module_name}/{self.job_id} is complete')
                break
            time.sleep(sleeping_time) #if job is not complete - wait some more
        job_report = f"*o{self.job_id.split('.')[0]}" #job stdout/stderr file name (job report)
        try:
            move(job_report, f'{self.output_path}/{self.module_name}_{self.job_name}_{self.job_id}.txt') #move job report to the output folder, where the rest of related files are generated
        except:
            print(f'Failed to move {job_report} to {self.output_path}')
        #os.system(f"mv {job_report} {self.output_path}/{self.module_name}_{self.job_name}_{self.job_id}.txt") #move job report to the output folder, where the rest of related files are generated


    def run_module_cluster(self, job_count): #AKA do_not_mess_with_my_quatation_marks
        '''
        Runs module on the login node of the HPC cluster, given number of jobs to run in parallel (int).
        Allows the snakemake to do job submissions to the computing nodes automatically.     
        '''
        #job_submission command to be used by snakmake to automatically submit jobs to HPC; stuff in curly brackets are snakemake arguments, not python variables
        if os.path.basename(self.cluster_config_path) == 'cluster.yaml':
            job_submission_command = '"qsub -N {cluster.jobname} -l nodes={cluster.nodes}:ppn={cluster.procs},pmem={cluster.pmem},walltime={cluster.walltime},feature={cluster.feature},file={cluster.file} -q {cluster.queue} -j {cluster.jobout} -o {cluster.outdir} -A {cluster.account} -V"'
        elif os.path.basename(self.cluster_config_path) == 'cluster_slurm.yaml':
            job_submission_command = '"sbatch --job-name {cluster.jobname} -N {cluster.nodes} --ntasks={cluster.ppn} --mem-per-cpu={cluster.mempc} -t {cluster.time} -o {cluster.outdir}{cluster.output} -e {cluster.outdir}{cluster.error} --export=ALL"'
        else:
            job_submission_command = '"qsub -N {cluster.jobname} -l procs={cluster.procs},pmem={cluster.pmem},walltime={cluster.walltime},feature={cluster.feature} -q {cluster.queue} -j {cluster.jobout} -o {cluster.outdir} -A {cluster.account} -V"'
        #shell command run by the wrapper (includes qsub command as substring);
        shell_command = f'''
        eval "$(conda shell.bash hook)";
        source activate /mnt/home/$(whoami)/.conda/envs/mamba_env/envs/snakemake;
        snakemake --reason --nolock --restart-times {self.retry_times} --jobs {job_count} --cluster-config {self.cluster_config_path} --cluster-status {self.status_script} --cluster-cancel qdel --configfile {self.config_file_path} --snakefile {self.snakefile_path} --keep-going --use-envmodules --use-conda --conda-frontend conda --rerun-incomplete --latency-wait 30 {self.force_all} {self.force_specific} {self.dry_run} --cluster {job_submission_command} {self.rule_graph} '''
        try:
            process_data = subprocess.check_call(shell_command, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as smk_error:
            smk_log            = smk_error.output
            failed_samples_tag = 'Out of jobs ready to be started, but not all files built yet.'

            if re.search(failed_samples_tag, smk_log):
                raise Exception(f"{self.module_name} module process: {failed_samples_tag}")
                #case 1: snakemake throws an error if it is out of jobs - workflow restart required
            else:
                #case 2: snakemake throws an error if there is a bug in the workflow code - fix required
                raise Exception(f"{self.module_name} module process running error: {smk_error.output}")
        except KeyboardInterrupt as ki:
            raise Exception(f"{self.module_name} was interrupted by the user: {ki}")
            #case 3 - keyboard interrupt by the user
        else:
            #if the workflow finished normally
            print(process_data)


    def pack_failed(self):
        '''
        Parses check_note_{self.module_name} of {self.sample_sheet} to get list of samples with at least 1 missing file.
        Creates {self.output_path}_failed_{self.module_name}_{timestamp} folder under parent folder of {self.output_path}.
        Moves all folders and files related to failed samples 
        from {self.output_path} to {self.output_path}_failed_{self.module_name}_{timestamp}.
        Sets {self.failed_stamp} to {timestamp} (default = None).
        '''
        failed_samples  = self.sample_sheet[self.sample_sheet[f'check_note_{self.module_name}'].str.contains('False')]['sample_id'].tolist()
        timestamp       = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
        failed_dir_path = f'{os.path.abspath(self.output_path)}_failed_{self.module_name}_{timestamp}'

        os.makedirs(failed_dir_path, mode=775, exist_ok=True)

        for id in failed_samples:
            for f in glob.glob(f'{self.output_path}{id}*'):
                move(f, f'{failed_dir_path}/')

        self.failed_stamp = timestamp


    def run_module(self, job_count, jobscript_path='./subscripts/ardetype_jobscript.sh'):
        '''Runs module on hpc as job or as snakemake submitter (on login node), based on self.run_mode value (True - job, False - submitter).'''
        if self.run_mode:
            self.submit_module_job(jobscript_path)
            self.check_job_completion()
        else:
            self.run_module_cluster(job_count)

        
    def add_taxonomy_column(self):
        '''Reads taxonomy information from self.aggr_taxonomy_path into self.taxonomy_dict 
        and adds taxonomy information as new column to the self.sample_sheet.'''
        self.taxonomy_dict = hk.read_json_dict(self.aggr_taxonomy_path)
        self.sample_sheet = hk.map_new_column(self.sample_sheet,self.taxonomy_dict,'sample_id','taxonomy')


    def clear_working_directory(self):
        '''Moves all files from working directory to source directory stored in self.cleanup_dict.'''
        # with open(f'{self.config_file["output_directory"]}/cleanup_dict.json', 'w+') as f:
        #     json.dump(self.cleanup_dict, f, indent=4)
        for key in self.cleanup_dict: 
            try:
                move(key, self.cleanup_dict[key])
            except:
                continue
            

    def files_to_wd(self, redirect_filter:dict=None):
        '''
        Moves all input files from input and output directories to working directory before running snakemake.
        If redirect_filter dictionary is passed, checks each files against the keys of it. If a match is found,
        sets source path in self.cleanup_dict to the value mapped to the corresponding key of redirect_filter (only one filter applied to each file).
        '''
        os.makedirs(os.path.abspath(self.config_file['work_dir']), exist_ok=True)
        for format in self.input_dict:
            map_dict = {}
            for source_path in self.input_dict[format]:
                full_path = f"{self.config_file['work_dir']}/{os.path.basename(source_path)}" #full input path to file
                if redirect_filter is not None: #if redirection was requested
                    for filter in redirect_filter: #starting to check filters against file names
                        if filter in source_path: #if match
                            map_dict[full_path] = redirect_filter[filter] #redirect
                            try:
                                # copy(source_path, os.path.abspath(self.config_file['work_dir']))
                                move(source_path, os.path.abspath(self.config_file['work_dir']))  #move to wd
                            except:
                                break #match found by moving file was not succesful
                            break #stop matching filters
                    if full_path not in map_dict: #if all filters are parsed but no match (if match happend, the full path will be in map_dict)
                        map_dict[full_path] = source_path #no redirection
                        try:
                            # copy(source_path, os.path.abspath(self.config_file['work_dir']))
                            move(source_path, os.path.abspath(self.config_file['work_dir']))
                        except Exception as e:
                            print(e)
                            continue
                else: #if no filtering required during function call
                    map_dict[full_path] = source_path
                    try:
                        move(source_path, os.path.abspath(self.config_file['work_dir']))
                    except:
                        continue
            self.cleanup_dict.update(map_dict) #add new entries to self.cleanup_dict - these are use to place files back to source/redirect location during cleanup


    def fold_output(self):
        '''Creates a folder for each sample_id in self.sample_sheet and self.removed_samples.
        Structures the pipeline output by putting all targets for each sample into curresponding folder.'''
        full_sample_list = self.sample_sheet['sample_id'].tolist() 
        if not self.removed_samples.empty: full_sample_list += self.removed_samples['sample_id'].to_list()
        for sample_id in full_sample_list: 
            os.makedirs(f'{self.output_path}folded_{sample_id}_output', exist_ok=True)
            for file in glob.glob(f'{self.output_path}{sample_id}*'): 
                try:
                    move(file, f'{self.output_path}folded_{sample_id}_output/')
                except Exception as e:
                    print(file, e)
                    continue

        # moving reports to separate directory - list of expected reports is defined in config_files/json/module_data/
        os.makedirs(f'{self.output_path}reports', exist_ok=True)
        outdir_listing = os.listdir(self.output_path)
        for suffix in module_data['reports']:
            for file in outdir_listing:
                if file.endswith(suffix):
                    try:
                        move(f'{self.output_path}{file}', f'{self.output_path}reports/')
                    except:
                        continue
        

    def unfold_output(self):
        '''Moves target files outside of folders created by fold_output method in order to avoid having to move file out manually to do a rerun.'''
        for file in glob.glob(f'{self.output_path}folded_*_output/*'):
            try:
                move(file, self.output_path)
            except:
                continue

        for file in glob.glob(f'{self.output_path}reports/*'):
            try:
                move(file, self.output_path)
            except:
                continue

    def set_permissions(self, permissions:str='775'):
        '''Given Linux permission string in numeric format, sets requested permissions (775 by default) recursively on the contents of self.output_path.'''
        hk.asign_perm_rec(self.output_path, permissions)

