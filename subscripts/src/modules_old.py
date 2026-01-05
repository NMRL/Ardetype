from .utilities import Housekeeper as hk
import os, warnings, re, subprocess, shutil, time, pandas as pd, glob, sys
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
            snakemake_cpus      : int,
            force_all           : bool, 
            pack_output         : bool, 
            rules_to_rerun      : list,
            run_local           : bool,
            nanopore_mode       : bool,
            fasta_mode          : bool
            ) -> None:
        
        self.run_mode            = run_mode #If true, snakemake will be run as single job, else - will run as job submitter on the login node
        self.run_local           = run_local #If true, each module will be executed on the machine where the wrapper is run
        self.job_id              = None  #Will be added if self.run_mode is True and job was submitted to HPC; filled by submit_module_job
        self.taxonomy_dict       = None   #Required if module creates different targets for different samples based on taxonomy information; filled by add_taxonomy_column
        self.module_name         = module_name #To be used in configuration file & sample_sheet file + to connect between modules (using remove_invalid_samples)
        self.input_path          = input_path #to the folder containing fasta/fastq.gz files
        self.output_path         = f"{os.path.abspath(output_path)}/" #Path to the output folder, where files will be saved (converted to full path)
        self.target_list         = None #List of all target files the module expects to create; filled by fill_target_list
        self.sample_sheet        = None #to store current state of sample_sheet dataframe; filled by create_sample_sheet; altered by fill_sample_sheet & receive_sample_sheet
        self.aggr_taxonomy_path  = f'{os.path.abspath(self.output_path)}/core_aggregated_taxonomy.json' #where to look for top kraken2 hits if snakemake will produce it; used by add_taxonomy_column
        self.config_file_path    = f'{os.path.abspath(self.output_path)}/config.yaml' #where to look for operational copy of the configuration file; used by submit_module_job & run_module_cluster
        self.cluster_config_path = cluster_config_path #where to look for job resource definition file; used by run_module_cluster
        self.config_file         = hk.read_yaml(module_config) if isinstance(module_config, str) else module_config #read module configuration from file if string is supplied (path expected); else - reads dictionary; used by add_module_targets, add_output_dir, write_module_config
        self.input_dict          = {} #to store input file paths for each file extension; used by fill_input_dict, fill_sample_sheet, add_fasta_samples
        self.patterns            = patterns #to store file extension patterns of expected input files; used by fill_input_dict; fill_sample_sheet
        self.job_name            = job_name #to store job name if self.run_mode is True; used by check_job_completion
        self.targets             = targets #to store file extensions of expected output files; used by fill_target_list, check_module_output
        self.requests            = requests #to store file extensions for files that are neccessary to run the modules; used by remove_invalid_samples
        self.snakefile_path      = snakefile_path #to the rule file to be run as single job on HPC if self.run_mode is True; used by submit_module_job
        self.retry_times         = retry_times #number of times snakemake will attempt to rerun failed jobs (default=3); used by run_module_cluster
        self.force_all           = "--forceall" if force_all else "" #to store forceall flag if it is supplied, else empty string is stored
        self.snakemake_cpus      = snakemake_cpus
        self.removed_samples     = pd.DataFrame() #to store dataframe containing information about samples that were deemed invalid by the module
        self.pack_output         = pack_output #switch to control putting output files into one folder named after sample_id; used by fold_output
        self.status_script       = f"{os.path.dirname(Path(__file__).parents[0].absolute())}/pbs-status.py"
        self.failed_stamp        = None #added if module has failed to produce requested files for 1 or more steps of the workflow
        self.rules_to_rerun      = [rule for rule in rules_to_rerun if rule in self._get_rule_names_from_snakefile(self.snakefile_path)] if rules_to_rerun is not None else []
        self.force_specific      = f"-R {' '.join(self.rules_to_rerun)}" if self.rules_to_rerun else ""
        self.fasta_mode          = fasta_mode
        self.nanopore_mode       = nanopore_mode
        
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


    def fill_input_dict(self, substring_list=['reads_unclassified', 'reads_classified'], mixed:bool=False, empty:bool=False, pattern_path:str=""):
        '''Fills self.input_dict using self.input_path and self.module_name by
        mapping each file format to the list of files of that format, found in the self.input_path, 
        excluding files that contain substrings in their names (supply None to avoid excluding files).
        If some files of required format are missing, raises an exception, indicating missing file format.'''
        if pattern_path:
            for fmt in self.patterns['inputs']:
                if fmt in pattern_path:
                    if isinstance(self.patterns['inputs'][fmt], list):
                        for pattern in self.patterns['inputs'][fmt]:
                            parsed_files = hk.parse_folder(
                                self.input_path,substr_lst=substring_list,
                                file_fmt_str=pattern)
                            if not parsed_files:
                                print(f'Missing {fmt} files in input directory', file=sys.stderr)
                                sys.exit(1)
                            self.input_dict[pattern] = parsed_files
                    elif isinstance(self.patterns['inputs'][fmt], str):
                        parsed_files = hk.parse_folder(
                                self.input_path,substr_lst=substring_list,
                                file_fmt_str=self.patterns['inputs'][fmt])
                        if not parsed_files:
                            print(f'Missing {fmt} files in input directory', file=sys.stderr)
                            sys.exit(1)
                        self.input_dict[self.patterns['inputs'][fmt]] = parsed_files

        elif mixed or empty:
            for fmt in self.patterns['inputs']['required']:
                self.input_dict[fmt] = hk.parse_folder(self.input_path,substr_lst=substring_list, file_fmt_str=fmt)
                if not self.input_dict[fmt]:
                    print(f'Missing {fmt} files in input directory', file=sys.stderr)
                    sys.exit(1)
            
            if not empty:
                for fmt in self.patterns['inputs']['optional']:
                    parsed_files = hk.parse_folder(self.input_path,substr_lst=substring_list, file_fmt_str=fmt)
                    if parsed_files:
                        self.input_dict[fmt] = parsed_files

    def fill_sample_sheet(self, pattern:str="001.fastq.gz"):
        '''
        Initializes self.sample_sheet to pandas dataframe, using self.input_dict and self.module_name (restricted to fastq & fasta inputs).
        '''
        ###Development - create sample sheet from paired or unpaired files; update sample sheet with paired or unpaired files
        if len(self.input_dict) < 2: #only one file extension is used - assumed fastq.gz
            self.sample_sheet = hk.create_sample_sheet(self.input_dict[pattern],self.patterns['sample_sheet'],mode=0)
        else: #fastq & fasta assumed
            self.sample_sheet = hk.create_sample_sheet(self.input_dict[pattern],self.patterns['sample_sheet'],mode=0)
            fasta_dict = {re.sub("_contigs.fasta","",os.path.basename(contig)):contig for contig in self.input_dict["_contigs.fasta"]}
            self.sample_sheet = hk.map_new_column(self.sample_sheet,fasta_dict,'sample_id','fa')


    def get_sample_groups(self, regexp:str=f'(_R[1,2]_001.fastq.gz|_ONT.fastq.gz)', fasta:bool=False):
        '''Extracts grouping information for samples based on the available input files'''
        if fasta:
            file_list = glob.glob(os.path.join(self.input_path, '*'))
            file_map = {}
            for f in file_list:
                sample_id = os.path.basename(re.sub(r'\_contigs.fasta', '', f))
                finding = re.search(r'\_contigs.fasta', f)
                if finding is not None:
                    finding = finding.group(0)
                    file_map[sample_id] = file_map.get(sample_id, list())
                    file_map[sample_id].append(finding)

            df = pd.DataFrame.from_dict(file_map, orient='index').astype(str).reset_index()
            df.columns = ['sample_id', 'FA']
            fa_case = ~df.FA.str.contains('None')

            ont_fa = df[fa_case]

            #Add groups to the sample sheet
            self.sample_sheet.loc[self.sample_sheet['sample_id'].isin(ont_fa['sample_id']), 'sample_group'] = 'FA'

            #Add mark file accordingly
            if not ont_fa.empty:
                f = open(os.path.join(self.output_path,'FA_mark'), 'w')
                f.close()
        else:
            file_list = glob.glob(os.path.join(self.input_path, '*'))
            file_map = {}
            for f in file_list:
                sample_id = os.path.basename(re.sub(regexp, '', f))
                finding = re.search(regexp, f)
                if finding is not None:
                    finding = finding.group(0)
                    file_map[sample_id] = file_map.get(sample_id, list())
                    file_map[sample_id].append(finding)
            #Equalize column length
            for k, v in file_map.items():
                ln = len(v)
                if ln < 3:
                    for _ in range(3 - ln):
                        file_map[k].append(str(None))
                file_map[k].sort()
            
            #Define groups of samples
            df = pd.DataFrame.from_dict(file_map, orient='index').reset_index()
            df.columns = ['sample_id', 'ONT', 'ILL1', 'ILL2']
            ont_case = df.ONT.str.contains('None') & df.ILL1.str.contains('None') & ~df.ILL2.str.contains('None')
            ill_case = df.ONT.str.contains('None') & ~df.ILL1.str.contains('None') & ~df.ILL2.str.contains('None')
            ful_case = ~df.ONT.str.contains('None') & ~df.ILL1.str.contains('None') & ~df.ILL2.str.contains('None')


            ont_pg = df[ont_case]
            ill_gp = df[ill_case]
            ful_gp = df[ful_case]


            #Add groups to the sample sheet
            self.sample_sheet.loc[self.sample_sheet['sample_id'].isin(ont_pg['sample_id']), 'sample_group'] = 'ONT'
            self.sample_sheet.loc[self.sample_sheet['sample_id'].isin(ill_gp['sample_id']), 'sample_group'] = 'ILL'
            self.sample_sheet.loc[self.sample_sheet['sample_id'].isin(ful_gp['sample_id']), 'sample_group'] = 'FUL'

            #Add mark file accordingly
            if not ont_pg.empty:
                f = open(os.path.join(self.output_path,'ONT_mark'), 'w')
                f.close()
            elif not ful_gp.empty:
                f = open(os.path.join(self.output_path,'HYB_IO_mark'), 'w')
                f.close()
            else:
                f = open(os.path.join(self.output_path,'ILL_mark'), 'w')
                f.close()


    def fill_target_list(self, taxonomy_based:bool=False, mixed:bool=False, empty:bool=False, grouped:bool=False):
        '''Fills self.target_list using data stored in self.sample_sheet instance variable.'''
        if taxonomy_based:#specific targets for each species
            self.target_list = [os.path.join(self.output_path,f'{id}{tmpl}') for idx, id in enumerate(self.sample_sheet['sample_id']) for tmpl in self.targets[self.sample_sheet['taxonomy'][idx]]]
        elif mixed or empty:#both species-specific and non-specific targets
            self.target_list = [os.path.join(self.output_path,f'{id}{tmpl}') for id in self.sample_sheet['sample_id'].to_list()+self.removed_samples['sample_id'].to_list() for tmpl in self.targets['general']]
            if not empty:
                self.target_list += [os.path.join(self.output_path,f'{id}{tmpl}') for idx, id in enumerate(self.sample_sheet['sample_id']) for tmpl in self.targets[self.sample_sheet['taxonomy'][idx]]]
        else:#only non-specific targets
            if grouped:
                self.target_list = [os.path.join(self.output_path,f'{id}{tmpl}') for idx, id in enumerate(self.sample_sheet['sample_id']) for tmpl in self.targets[self.sample_sheet['sample_group'][idx]]]
            else:
                self.target_list = [os.path.join(self.output_path,f'{id}{tmpl}') for id in self.sample_sheet['sample_id'] for tmpl in self.targets]
            

    def make_output_dir(self):
        '''Creates output directory (if not present in the file system) using self.output_path.'''
        os.makedirs(self.output_path, exist_ok=True)


    def write_sample_sheet(self):
        '''Creates sample_sheet.csv file in the self.output_path folder, using self.sample_sheet.'''
        self.sample_sheet.to_csv(os.path.join(self.output_path, "sample_sheet.csv"), header=True, index=False)


    def add_module_targets(self):
        '''Updates self.config_file, using self.module_name.'''
        validation_code = hk.validate_yaml(self.config_file)
        if validation_code != 0:
            print(f'Config file structure differs from expected template (see config_files/yaml/config_modular_local.yaml for example).', file=sys.stderr)
            sys.exit(1)

        output_code = hk.edit_nested_dict(config_dict=self.config_file, param=f"{self.module_name}_target_files", new_value=self.target_list)
        if output_code != 0:
            nl = '\n'
            print(f'Failed to set {self.module_name}_target_files value to\n\n{nl.join(self.target_list)}\n\nPlease ensure that config file contains {self.module_name}_target_files parameter (see config_files/yaml/config_modular_local.yaml for example).', file=sys.stderr)
            sys.exit(1)


    def add_output_dir(self):
        '''Updates self.config_file using self.output_path.'''

        validation_code = hk.validate_yaml(self.config_file)
        if validation_code != 0:
            print(f'Config file structure differs from expected template (see config_files/yaml/config_modular_local.yaml for example).', file=sys.stderr)
            sys.exit(1)

        output_code = hk.edit_nested_dict(config_dict=self.config_file, param="output_directory", new_value=self.output_path)
        if output_code != 0:
            print(f'Failed to set output_directory value to {self.output_path}.\nPlease ensure that config file contains output_directory parameter (see config_files/yaml/config_modular_local.yaml for example).', file=sys.stderr)
            sys.exit(1)

        inpath = self.output_path if self.output_path == self.input_path else os.path.abspath(self.input_path) + "/"
        output_code = hk.edit_nested_dict(config_dict=self.config_file, param="input_directory", new_value=inpath)
        if output_code != 0:
            print(f'Failed to set input_directory value to {self.input_path}.\nPlease ensure that config file contains input_directory parameter (see config_files/yaml/config_modular_local.yaml for example).', file=sys.stderr)
            sys.exit(1)

        #Updating paths in config
        ardetype_path = str(Path(os.path.abspath('./')))
        for p in ["databases"]:
            if not os.path.isabs(self.config_file[p]):
                self.config_file[p] = os.path.join(ardetype_path, self.config_file[p])

        for p in ["scratch"]:
            if not os.path.isabs(self.config_file[p]):
                self.config_file[p] = os.path.join(ardetype_path, self.output_path, self.config_file[p])

        for substr in ["_database", "_db", "_sif"]:
            hk.join_sif_paths(self.config_file, substr, ardetype_path)


    def write_module_config(self):
        '''Writes self.config_file to the self.output_path'''
        hk.write_yaml(self.config_file, os.path.join(self.output_path,'config.yaml'))


    def check_module_output(self, mixed:bool=False):
        '''Checks if output files are generated according to self.module_name and adds check_note_{self.module_name} column 
        to the self.sample_sheet dataframe, where boolean value is stored for each expected file.'''
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
    


    def run_modules_local(self):
        input_dir = Path(self.output_path)
        failed_dir = Path(f"{os.path.dirname(self.output_path)}_failed")
        failed_samples = []

        def remaining_samples(class_instance):
            if class_instance.fasta_mode:
                return {p.name.split("_contigs.fasta")[0] for p in input_dir.glob("*_contigs.fasta")}
            elif class_instance.nanopore_mode:
                return {p.name.split("_ONT.fastq.gz")[0] for p in input_dir.glob("*_ONT.fastq.gz")}
            else:
                return {p.name.split("_R1_001.fastq.gz")[0] for p in input_dir.glob("*_R1_001.fastq.gz")}

        while True:
            try:
                if failed_samples:
                    if self.nanopore_mode:
                        self.fill_input_dict(pattern_path='ONT')
                        self.fill_sample_sheet(pattern=self.patterns['inputs']['ONT'])
                    elif self.fasta_mode:
                        self.fill_input_dict(pattern_path='FASTA')
                        self.fill_sample_sheet(pattern=self.patterns['inputs']['FASTA'])
                    else:
                        self.fill_input_dict(pattern_path='ILL')
                        self.fill_sample_sheet(pattern=self.patterns['inputs']['ILL'])
                    self.get_sample_groups(fasta=self.fasta_mode)
                    self.write_sample_sheet()
                    self.fill_target_list(grouped=True)
                    self.add_module_targets()
                    self.add_output_dir()
                    self.write_module_config()
                cmd = f'''
                source "$(conda info --base)/etc/profile.d/conda.sh"
                conda activate $(dirname $(dirname $(which conda)))/envs/ardetype
                snakemake --cores {self.snakemake_cpus} --reason --nolock --restart-times {self.retry_times} --resources API_calls=1 --configfile {self.config_file_path} --snakefile {self.snakefile_path} --keep-going --rerun-incomplete --latency-wait 30 {self.force_all} {self.force_specific} | tee /tmp/snakemake_run_$(date +%s).log
                '''
                result = subprocess.run(cmd, shell=True, executable="/bin/bash", check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                self.job_id = result.stdout
                break

            except subprocess.CalledProcessError as e:
                failed_samples = hk.get_failed_sample_counts(e.stderr)
                failed_dir.mkdir(exist_ok=True)

                if not failed_samples:
                    print(f"Snakemake failed, but no failed samples found: {e.stderr}", file=sys.stderr)
                    sys.exit(1)

                for sample in failed_samples:
                    if failed_samples[sample] > int(self.retry_times):
                        for f in input_dir.glob(f"{sample}*"):
                            f.rename(failed_dir / f.name)

                if not remaining_samples(class_instance=self):
                    print(f"Missing input files in {input_dir} after moving failed ones to {failed_dir}, exiting.", file=sys.stderr)
                    sys.exit(1)
            except KeyboardInterrupt:
                print("Interrupted by user.", file=sys.stderr)
                sys.exit(1)
            finally:
                self.clear_scratch()

    # def run_modules_local(self):
    #     try:
    #         cmd = f'''
    #         source "$(conda info --base)/etc/profile.d/conda.sh"
    #         conda activate $(dirname $(dirname $(which conda)))/envs/ardetype
    #         snakemake --cores {self.snakemake_cpus} --reason --nolock --restart-times {self.retry_times} --resources API_calls=1 --configfile {self.config_file_path} --snakefile {self.snakefile_path} --keep-going --rerun-incomplete --latency-wait 30 {self.force_all} {self.force_specific}
    #         '''
    #         result = subprocess.run(cmd, shell=True, executable="/bin/bash", check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    #         self.job_id = result.stdout
    #     except KeyboardInterrupt:
    #         print(f"{self.module_name} interrupted by SIGTERM : Run the same command to continue (make sure to include timestamp added to batch folder name).", file=sys.stderr)
    #         sys.exit(1)
    #     except Exception as e:
    #         print(f"{self.module_name} finished with runtime error:\nSTDOUT:\n{e.stdout}\nSTDERR:\n", file=sys.stderr)
    #         sys.exit(1)
    #     finally:
    #         self.clear_scratch()


    def submit_module_job(self, jobscript_path):
        """
        Submits module as a job to HPC cluster, given path to the job_script, setting self.job_id to the bytestring representing job_id.
        All rules are run sequentially using same set of resources.
        """
        shutil.copy(jobscript_path, f'{self.output_path}ardetype_jobscript.sh') #to avoid running the template file

        try:
            self.job_id = subprocess.check_output(['qsub', '-F', f'{self.snakefile_path} {self.config_file_path}', f'{self.output_path}ardetype_jobscript.sh'])
            os.remove(f'{self.output_path}ardetype_jobscript.sh') #cleanup

        except subprocess.CalledProcessError as msg:
            print(f"{self.module_name} module job submission failed with error:{msg}", file=sys.stderr)
            sys.exit(1)
        finally:
            self.clear_scratch()
        

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
        cluster = hk.read_yaml(self.cluster_config_path)['cluster']
        self.job_cancel = 'qdel'
        if cluster == 'pbs':
            job_submission_command = '"qsub -N {cluster.jobname} -l nodes={cluster.nodes}:ppn={cluster.procs},pmem={cluster.pmem},walltime={cluster.walltime},feature={cluster.feature},file={cluster.file} -q {cluster.queue} -j {cluster.jobout} -o {cluster.outdir} -A {cluster.account} -V"'
        elif cluster == 'slurm':
            job_submission_command = '"sbatch --parsable --job-name {cluster.jobname} -N {cluster.nodes} --ntasks={cluster.npn} --mem-per-cpu={cluster.mempc} -t {cluster.time} -o {cluster.outdir}{cluster.output} -e {cluster.outdir}{cluster.error} --export=ALL"'
            self.status_script = f"{os.path.dirname(Path(__file__).parents[0].absolute())}/slurm-status.py"
            self.job_cancel = 'scancel'
        else:
            job_submission_command = '"qsub -N {cluster.jobname} -l procs={cluster.procs},pmem={cluster.pmem},walltime={cluster.walltime},feature={cluster.feature} -q {cluster.queue} -j {cluster.jobout} -o {cluster.outdir} -A {cluster.account} -V"'
        #shell command run by the wrapper (includes qsub command as substring);
                # eval "$(conda shell.bash hook)";
        # source activate /mnt/home/$(whoami)/.conda/envs/mamba_env/envs/snakemake;
        shell_command = f'''
        snakemake --scheduler greedy --reason --nolock --restart-times {self.retry_times} --resources API_calls=1 --jobs {job_count} --cluster-config {self.cluster_config_path} --cluster-status {self.status_script} --cluster-cancel {self.job_cancel} --configfile {self.config_file_path} --snakefile {self.snakefile_path} --keep-going --use-envmodules --use-conda --conda-frontend conda --rerun-incomplete --latency-wait 300 {self.force_all} {self.force_specific} --cluster {job_submission_command}'''
        try:
            process_data = subprocess.check_call(shell_command, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as smk_error:
            smk_log            = smk_error.output
            failed_samples_tag = 'Out of jobs ready to be started, but not all files built yet.'

            if re.search(failed_samples_tag, smk_log):
                #case 1: snakemake throws an error if it is out of jobs - workflow restart required
                print(f"{self.module_name} module failed for {failed_samples_tag} samples.", file=sys.stderr)
                sys.exit(1)
            else:
                #case 2: snakemake throws an error if there is a bug in the workflow code - fix required
                print(f"{self.module_name} module code contains error: {smk_error.output}", file=sys.stderr)
                sys.exit(1)
        except KeyboardInterrupt as ki:
            #case 3 - keyboard interrupt by the user
            print(f"{self.module_name} was interrupted by the user: {ki}", file=sys.stderr)
            sys.exit(1)           
        else:
            #if the workflow finished normally
            print(process_data)


    def run_module(self, job_count, jobscript_path='./subscripts/ardetype_jobscript.sh'):
        '''Runs module on hpc as job or as snakemake submitter (on login node), based on self.run_mode value (True - job, False - submitter).'''
        if not self.run_local:
            if self.run_mode:
                self.submit_module_job(jobscript_path)
                self.check_job_completion()
            else:
                self.run_module_cluster(job_count)
        else:
            self.run_modules_local()
        
    def add_taxonomy_column(self):
        '''Reads taxonomy information from self.aggr_taxonomy_path into self.taxonomy_dict 
        and adds taxonomy information as new column to the self.sample_sheet.'''
        if 'taxonomy' in self.sample_sheet.columns:
            self.sample_sheet = self.sample_sheet.drop('taxonomy', axis=1)
        self.taxonomy_dict = hk.read_json_dict(self.aggr_taxonomy_path)
        self.sample_sheet = hk.map_new_column(self.sample_sheet,self.taxonomy_dict,'sample_id','taxonomy')
        
    def clear_scratch(self):
        ''' Removed every file and folder from scratch directory.'''
        if os.path.isdir(self.config_file["scratch"]):
            dir_contents = [os.path.join(self.config_file["scratch"],p) for p in os.listdir(self.config_file["scratch"])]
            for path in dir_contents:
                if os.path.isdir(path):
                    shutil.rmtree(path)
                else:
                    os.remove(path)


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
        unfolding_path = self.output_path if self.input_path == self.output_path else self.input_path
        
        for file in glob.glob(f'{unfolding_path}folded_*_output/*'):
            try:
                move(file, unfolding_path)
            except Exception as e:
                print(e)
                continue

        for file in glob.glob(f'{unfolding_path}reports/*'):
            try:
                move(file, unfolding_path)
            except Exception as e:
                print(e)
                continue

    def set_permissions(self, permissions:str='775'):
        '''Given Linux permission string in numeric format, sets requested permissions (775 by default) recursively on the contents of self.output_path.'''
        hk.asign_perm_rec(self.output_path, permissions)

