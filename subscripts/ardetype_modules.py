from ardetype_utilities import *
import os, warnings, re, sys, subprocess, shutil, time
from itertools import chain
from getpass import getuser
from pathlib import Path
sys.path.insert(0, os.path.dirname(Path(__file__).absolute()))

#Suppressing pandas warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)


#Reading data used to build module objects
module_data = read_json_dict(f'{os.path.dirname(Path(__file__).parents[0].absolute())}/config_files/json/module_data.json')


####################
# Class definition
####################

class Module:
    '''Class represents single module of the ardetype pipeline'''

    def __init__(self, module_name: str, input_path: str, module_config, output_path: str, run_mode: bool, job_name: str, patterns: dict, targets: list, requests: dict, snakefile_path: str, cluster_config_path: str, dry_run: bool, force_all: bool, rule_graph: bool, pack_output:bool, unpack_output:bool) -> None:
        self.run_mode = run_mode #If true, snakemake will be run as single job, else - will run as job submitter on the login node
        self.job_id = None  #Will be added if self.run_mode is True and job was submitted to HPC; filled by submit_module_job
        self.taxonomy_dict = None   #Required if module creates different targets for different samples based on taxonomy information; filled by add_taxonomy_column
        self.module_name = module_name #To be used in configuration file & sample_sheet file + to connect between modules (using remove_invalid_samples)
        self.input_path = input_path #to the folder containing fasta/fastq.gz files
        self.output_path = f"{os.path.abspath(output_path)}/" #Path to the output folder, where files will be saved (converted to full path)
        self.target_list = None #List of all target files the module expects to create; filled by fill_target_list
        self.sample_sheet = None #to store current state of sample_sheet dataframe; filled by create_sample_sheet; altered by fill_sample_sheet & receive_sample_sheet
        self.aggr_taxonomy_path = f'{os.path.abspath(self.output_path)}/{self.module_name}_aggregated_taxonomy.json' #where to look for top kraken2 hits if snakemake will produce it; used by add_taxonomy_column
        self.config_file_path = f'{os.path.abspath(self.output_path)}/config.yaml' #where to look for operational copy of the configuration file; used by submit_module_job & run_module_cluster
        self.cluster_config_path = cluster_config_path #where to look for job resource definition file; used by run_module_cluster
        self.config_file = read_yaml(module_config) if isinstance(module_config, str) else module_config #read module configuration from file if string is supplied (path expected); else - reads dictionary; used by add_module_targets, add_output_dir, write_module_config
        self.input_dict = {} #to store input file paths for each file extension; used by fill_input_dict, fill_sample_sheet, add_fasta_samples
        self.patterns = patterns #to store file extension patterns of expected input files; used by fill_input_dict; fill_sample_sheet
        self.job_name = job_name #to store job name if self.run_mode is True; used by check_job_completion
        self.targets = targets #to store file extensions of expected output files; used by fill_target_list, check_module_output
        self.requests = requests #to store file extensions for files that are neccessary to run the modules; used by remove_invalid_samples
        self.snakefile_path = snakefile_path #to the rule file to be run as single job on HPC if self.run_mode is True; used by submit_module_job
        self.dry_run = "-np" if dry_run else "" #to store dry-run flag if it is supplied, else empty string is stored
        self.force_all = "--forceall" if force_all else "" #to store forceall flag if it is supplied, else empty string is stored
        self.rule_graph = f"--rulegraph | dot -Tpdf > {self.module_name}.pdf" if rule_graph else "" #to store rule_graph flag if it is supplied, else empty string is stored
        self.unpack_output = unpack_output #used to move files outside sample folders and do a rerun; used by unfold_output
        self.removed_samples = pd.DataFrame() #to store dataframe containing information about samples that were deemed invalid by the module
        self.pack_output = pack_output #switch to control putting output files into one folder named after sample_id; used by fold_output
        self.cleanup_dict = {} #to map origin paths of input files to path in working directory; filled by move_to_wd; used by clear_working_directory


    def fill_input_dict(self, substring_list=['reads_unclassified', 'reads_classified']):
        '''Fills self.input_dict using self.input_path and self.module_name by
        mapping each file format to the list of files of that format, found in the self.input_path, 
        excluding files that contain substrings in their names (supply None to avoid excluding files).
        If some files of required format are missing, raises an exception, indicating missing file format.'''
        for format in self.patterns['inputs']: 
            self.input_dict[format] = parse_folder(self.input_path,substr_lst=substring_list, file_fmt_str=format)
            if not self.input_dict[format]: raise Exception(f'Missing {format} files in input directory')
   

    def create_sample_sheet(self):
        '''
        Initializes self.sample_sheet to pandas dataframe, using self.input_dict and self.module_name (restricted to fastq & fasta inputs).
        '''
        if len(self.input_dict) < 2: #only one file extension is used - assumed fastq.gz
            self.sample_sheet = create_sample_sheet(self.input_dict["001.fastq.gz"],self.patterns['sample_sheet'],mode=0)
        else: #fastq & fasta assumed
            self.sample_sheet = create_sample_sheet(self.input_dict["001.fastq.gz"],self.patterns['sample_sheet'],mode=0)
            fasta_dict = {re.sub("_contigs.fasta","",os.path.basename(contig)):contig for contig in self.input_dict["_contigs.fasta"]}
            self.sample_sheet = map_new_column(self.sample_sheet,fasta_dict,'sample_id','fa')


    def fill_target_list(self, taxonomy_based=False):
        '''Fills self.target_list using data stored in self.sample_sheet instance variable.'''
        if taxonomy_based:#specific targets for each species
            self.target_list = [f'{self.output_path}{id}{tmpl}' for idx, id in enumerate(self.sample_sheet['sample_id']) for tmpl in self.targets[self.sample_sheet['taxonomy'][idx]]]
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
        output_code = edit_nested_dict(config_dict=self.config_file, param=f"{self.module_name}_target_files", new_value=self.target_list)
        validation_code = validate_yaml(self.config_file)
        if not output_code == 0: raise Exception(f'Config editing failed with error code {output_code}')
        elif not validation_code == 0: raise Exception(f'Config validation failed with error code {validation_code}')


    def add_output_dir(self):
        '''Updates self.config_file using self.output_path.'''
        output_code = edit_nested_dict(config_dict=self.config_file, param="output_directory", new_value=self.output_path)
        validation_code = validate_yaml(self.config_file)
        if not output_code == 0: raise Exception(f'Config editing failed with error code {output_code}')
        elif not validation_code == 0: raise Exception(f'Config validation failed with error code {validation_code}')


    def write_module_config(self):
        '''Writes self.config_file to the self.output_path'''
        write_yaml(self.config_file, f'{self.output_path}config.yaml')


    def check_module_output(self):
        '''Checks if output files are generated according to self.module_name and adds check_note_{self.module_name} column 
        to the self.sample_sheet dataframe, where boolean value is stored for each expected file.'''
        ###Development - Automatically scale dirs_up depending on input structure - currently two dirs up max
        check_dict = check_file_existance(file_list=self.target_list)
        id_check_dict = {id:"" for id in self.sample_sheet['sample_id']}
        for file in check_dict:
            two_dirs_up = os.path.basename(os.path.dirname(os.path.dirname(file)))+"/"+os.path.basename(os.path.dirname(file))+"/"+os.path.basename(file) #required for outputs where directory patters are defined in addition to file extensions
            if isinstance(self.targets, list): #if only un-specific targets are supplied
                id = os.path.basename(re.sub("("+"|".join(self.targets)+")","",two_dirs_up))
            elif isinstance(self.targets, dict): #if taxonomy-based targets are supplied - all species-specific target lists are to be merged into one list using chain.from_iterables
                id = os.path.basename(re.sub("("+"|".join(chain.from_iterable(self.targets.values()))+")","",two_dirs_up))
            id_check_dict[id] += f"|{file}:{check_dict[file]}"
        self.sample_sheet = map_new_column(self.sample_sheet, id_check_dict, 'sample_id', f"check_note_{self.module_name}")


    def supply_sample_sheet(self): #getter, may not be required now as all variables are public, but makes it easier to encapsulate later, if needed
        '''Returns self.sample_sheet dataframe object.'''
        return self.sample_sheet


    def receive_sample_sheet(self, sample_sheet): #setter, may not be required now as all variables are public, but makes it easier to encapsulate later, if needed
        '''Inializes self.sample_sheet with external sample_sheet dataframe (used to connect modules).'''
        self.sample_sheet = sample_sheet


    def add_fasta_samples(self):
        '''Adds fa column with _contigs.fasta files to the self.sample_sheet dataframe.'''
        fasta_dict = {re.sub("_contigs.fasta","",os.path.basename(contig)):contig for contig in self.input_dict["_contigs.fasta"]}
        self.sample_sheet = map_new_column(self.sample_sheet,fasta_dict,'sample_id','fa')


    def remove_invalid_samples(self, connect_from_module_name, taxonomy_only=False):
        '''
        Removes samples that lack files, required by the current module, given supplier module name.
        If all samples are removed, returns 1 (int).
        '''
        if not taxonomy_only:
            if self.requests['check'] is not None: #if module requires certain file types to run rules that are not taxonomy-specific
                if isinstance(self.requests['check'],str): #if only one requirement is provided as string
                    self.removed_samples = self.sample_sheet[~self.sample_sheet[f'check_note_{connect_from_module_name}'].str.contains(self.requests['check'])].reset_index(drop=True)
                    self.sample_sheet = self.sample_sheet[self.sample_sheet[f'check_note_{connect_from_module_name}'].str.contains(self.requests['check'])].reset_index(drop=True)
                else:
                    for request in self.requests['check']: #if many requirements are provided as list of stings
                        self.removed_samples = self.sample_sheet[~self.sample_sheet[f'check_note_{connect_from_module_name}'].str.contains(request)].reset_index(drop=True)
                        self.sample_sheet = self.sample_sheet[self.sample_sheet[f'check_note_{connect_from_module_name}'].str.contains(request)].reset_index(drop=True)
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
        if not self.removed_samples.empty: self.removed_samples.to_csv(f"{self.output_path}removed_samples_{self.module_name}.csv", header=True, index=False)
            

    def submit_module_job(self, jobscript_path):
        """
        Submits module as a job to HPC cluster, given path to the job_script, setting self.job_id to the bytestring representing job_id.
        All rules are run sequentially using same set of resources.
        """
        shutil.copy(jobscript_path, f'{self.output_path}ardetype_jobscript.sh') #to avoid running the template file

        try:
            self.job_id = subprocess.check_output(['qsub', '-F', f'{self.snakefile_path} {self.config_file_path}', f'{self.output_path}ardetype_jobscript.sh'])
            os.system(f"rm {self.output_path}ardetype_jobscript.sh") #cleanup
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
        os.system(f"mv {job_report} {self.output_path}/{self.module_name}_{self.job_name}_{self.job_id}.txt") #move job report to the output folder, where the rest of related files are generated


    def run_module_cluster(self, job_count): #AKA do_not_mess_with_my_quatation_marks
        '''
        Runs module on the login node of the HPC cluster, given number of jobs to run in parallel (int).
        Allows the snakemake to do job submissions to the computing nodes automatically.     
        '''
        #job_submission command to be used by snakmake to automatically submit jobs to HPC; stuff in curly brackets are snakemake arguments, not python variables
        if os.path.basename(self.cluster_config_path) == 'cluster.yaml':
            job_submission_command = '"qsub -N {cluster.jobname} -l nodes={cluster.nodes}:ppn={cluster.ppn},pmem={cluster.pmem},walltime={cluster.walltime} -q {cluster.queue} -j {cluster.jobout} -o {cluster.outdir} -V"'
        elif os.path.basename(self.cluster_config_path) == 'cluster_slurm.yaml':
            job_submission_command = '"sbatch --job-name {cluster.jobname} -N {cluster.nodes} --ntasks={cluster.ppn} --mem-per-cpu={cluster.mempc} -t {cluster.time} -o {cluster.outdir}{cluster.output} -e {cluster.outdir}{cluster.error} --export=ALL"'
        #shell command run by the wrapper (includes qsub command as substring); to run in dry-run mode, add -np at the end of the snakemake command
        shell_command = f'''
        eval "$(conda shell.bash hook)";
        conda activate /mnt/home/$(whoami)/.conda/envs/mamba_env/envs/snakemake; 
        snakemake --jobs {job_count} --cluster-config {self.cluster_config_path} --cluster-cancel qdel --configfile {self.config_file_path} --snakefile {self.snakefile_path} --keep-going --use-envmodules --use-conda --conda-frontend conda --rerun-incomplete --latency-wait 30 {self.force_all} {self.rule_graph} {self.dry_run} --cluster {job_submission_command} '''
        print(shell_command)
        try:
            subprocess.check_call(shell_command, shell=True)
        except subprocess.CalledProcessError as msg:
            raise Exception(f"{self.module_name} module process running error: {msg}")


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
        self.taxonomy_dict = read_json_dict(self.aggr_taxonomy_path)
        self.sample_sheet = map_new_column(self.sample_sheet,self.taxonomy_dict,'sample_id','taxonomy')


    def clear_working_directory(self):
        """Moves all files from working directory to source directory stored in self.cleanup_dict."""
        for key in self.cleanup_dict: os.system(f"mv -n {key} {self.cleanup_dict[key]} 2> /dev/null")
            

    def files_to_wd(self):
        '''Moves all input files from input and output directories to working directory before running snakemake.'''
        os.makedirs(os.path.abspath(self.config_file['work_dir']), exist_ok=True)
        for format in self.input_dict:
            map_dict = {}
            for source_path in self.input_dict[format]:
                map_dict[f"{self.config_file['work_dir']}/{os.path.basename(source_path)}"] = source_path
                os.system(f"mv -n {source_path} {os.path.abspath(self.config_file['work_dir'])} 2> /dev/null")
            self.cleanup_dict.update(map_dict)


    def fold_output(self):
        '''Creates a folder for each sample_id in self.sample_sheet and self.removed_samples.
        Structures the pipeline output by putting all tartets for each sample into curresponding folder.'''
        full_sample_list = self.sample_sheet['sample_id'].tolist() 
        if not self.removed_samples.empty: full_sample_list += self.removed_samples['sample_id'].to_list()
        for sample_id in full_sample_list: 
            os.makedirs(f'{self.output_path}folded_{sample_id}_output', exist_ok=True)
            os.system(f'mv -n {self.output_path}{sample_id}* {self.output_path}folded_{sample_id}_output/ 2> /dev/null')


    def unfold_output(self):
        '''Moves target files outside of folders created by fold_output method in order to avoid having to move file out manually to do a rerun.'''
        for id in self.sample_sheet['sample_id']: os.system(f'mv -n {self.output_path}folded_{id}_output/* {self.output_path} 2> /dev/null')


###############################################
# Defining wrapper functions to call from main
###############################################

def run_all(args, num_jobs):
    '''Wrapper function to run all modules sequentially.'''
    core = Module(
            module_name='core',
            input_path=args.input,
            module_config=args.config,
            output_path=args.output_dir,
            run_mode=args.submit_modules,
            dry_run=args.dry_run,
            force_all=args.force_all,
            rule_graph=args.rule_graph,
            pack_output=args.pack_output,
            unpack_output=args.unpack_output,
            job_name=module_data['core']['job_name'],
            patterns=module_data['core']['patterns'],
            targets=module_data['core']['targets'],
            requests=module_data['core']['requests'],
            snakefile_path=module_data['snakefiles']['core'],
            cluster_config_path=module_data['cluster_config']
            )
    shell = Module(
        module_name='shell', 
        input_path=core.output_path, 
        module_config=core.config_file, 
        output_path=args.output_dir, 
        run_mode=args.submit_modules,
        dry_run=args.dry_run,
        force_all=args.force_all,
        rule_graph=args.rule_graph,
        pack_output=args.pack_output,
        unpack_output=args.unpack_output,
        job_name=module_data['shell']['job_name'],
        patterns=module_data['shell']['patterns'],
        targets=module_data['shell']['targets'],
        requests=module_data['shell']['requests'],
        snakefile_path=module_data['snakefiles']['shell'],
        cluster_config_path=module_data['cluster_config']
        )
    tip = Module(
        module_name='tip', 
        input_path=core.output_path,
        module_config=shell.config_file, 
        output_path=args.output_dir, 
        run_mode=args.submit_modules,
        dry_run=args.dry_run,
        force_all=args.force_all,
        rule_graph=args.rule_graph,
        pack_output=args.pack_output,
        unpack_output=args.unpack_output,
        job_name=module_data['tip']['job_name'],
        patterns=module_data['tip']['patterns'],
        targets=module_data['tip']['targets'],
        requests=module_data['tip']['requests'],
        snakefile_path=module_data['snakefiles']['tip'],
        cluster_config_path=module_data['cluster_config']
    )

    #Running core
    core.fill_input_dict()
    core.fill_sample_sheet()
    if core.unfold_output: core.unfold_output()
    core.make_output_dir()
    core.write_sample_sheet()
    core.fill_target_list()
    core.add_module_targets()
    core.add_output_dir()
    core.write_module_config()
    core.files_to_wd()
    try:
        core.run_module(job_count=num_jobs)
    except Exception as e:
        core.clear_working_directory() #to avoid manually moving files back to input
        raise e
    core.check_module_output()
    try:
        core.add_taxonomy_column()
    except FileNotFoundError as e: #it should be raised in dry-run mode as rule all of bact_core is not executed
        if core.dry_run == "" and core.rule_graph == "":
            raise e
    core.write_sample_sheet()
    core.clear_working_directory()

  
    #Connecting core to shell
    shell.receive_sample_sheet(core.supply_sample_sheet())
    if shell.dry_run == "" and shell.rule_graph == "": 
        samples_cleared = shell.remove_invalid_samples(connect_from_module_name='core') #in dry run mode none of the rules are executed, hence all samples will be removed, causing error
        shell.save_removed()
        if samples_cleared == 1: 
            if core.pack_output: core.fold_output()
            raise Exception('Missing files requested by bact_shell.')

    #Running shell
    shell.fill_input_dict()
    shell.add_fasta_samples()
    shell.write_sample_sheet()
    shell.fill_target_list()
    shell.add_module_targets()
    shell.write_module_config()
    shell.files_to_wd()
    try:
        shell.run_module(job_count=num_jobs)
    except Exception as e:
        shell.clear_working_directory()
        raise e
    shell.check_module_output()
    shell.write_sample_sheet()
    shell.clear_working_directory()

    # Connecting core to tip
    tip.receive_sample_sheet(shell.supply_sample_sheet())
    samples_cleared = tip.remove_invalid_samples(connect_from_module_name='core')
    tip.save_removed()
    if samples_cleared == 1: 
        if shell.pack_output: shell.fold_output()
        raise Exception('Missing files requested by bact_tip.')

    # Running tip
    tip.fill_input_dict(substring_list=None)
    tip.add_fasta_samples()
    tip.write_sample_sheet()
    tip.fill_target_list(taxonomy_based=True)
    tip.add_module_targets()
    tip.write_module_config()
    tip.files_to_wd()
    try:
        tip.run_module(job_count=num_jobs)
    except Exception as e:
        tip.clear_working_directory()
        raise e
    tip.check_module_output()
    tip.write_sample_sheet()
    tip.clear_working_directory()
    if tip.pack_output: tip.fold_output()


def run_core(args, num_jobs):
    '''Wrapper function to run only core module.'''
    core = Module(
        module_name='core',
        input_path=args.input,
        module_config=args.config,
        output_path=args.output_dir,
        run_mode=args.submit_modules,
        dry_run=args.dry_run,
        force_all=args.force_all,
        rule_graph=args.rule_graph,
        pack_output=args.pack_output,
        unpack_output=args.unpack_output,
        job_name=module_data['core']['job_name'],
        patterns=module_data['core']['patterns'],
        targets=module_data['core']['targets'],
        requests=module_data['core']['requests'],
        snakefile_path=module_data['snakefiles']['core'],
        cluster_config_path=module_data['cluster_config']
    )
    core.fill_input_dict()
    core.fill_sample_sheet()
    if core.unfold_output: core.unfold_output()
    core.make_output_dir()
    core.write_sample_sheet()
    core.fill_target_list()
    core.add_module_targets()
    core.add_output_dir()
    core.write_module_config()
    core.files_to_wd()
    try:
        core.run_module(job_count=num_jobs)
    except Exception as e:
        core.clear_working_directory()
        raise e
    core.check_module_output()
    try:
        core.add_taxonomy_column()
    except FileNotFoundError as error:
        if core.dry_run == "" and core.rule_graph == "":
            raise error
    core.write_sample_sheet()
    core.clear_working_directory()
    if core.pack_output: core.fold_output()


def run_shell(args, num_jobs):
    '''Wrapper function to run only shell module.'''
    shell = Module(
        module_name='shell', 
        input_path=args.input,
        module_config=args.config, 
        output_path=args.output_dir, 
        run_mode=args.submit_modules,
        dry_run=args.dry_run,
        force_all=args.force_all,
        rule_graph=args.rule_graph,
        pack_output=args.pack_output,
        unpack_output=args.unpack_output,
        job_name=module_data['shell']['job_name'],
        patterns=module_data['shell']['patterns'],
        targets=module_data['shell']['targets'],
        requests=module_data['shell']['requests'],
        snakefile_path=module_data['snakefiles']['shell'],
        cluster_config_path=module_data['cluster_config']
    )
    shell.fill_input_dict()
    shell.fill_sample_sheet()
    if shell.unpack_output: shell.unfold_output()
    shell.make_output_dir()
    shell.write_sample_sheet()
    shell.fill_target_list()
    shell.add_module_targets()
    shell.add_output_dir()
    shell.write_module_config()
    shell.files_to_wd()
    try:
        shell.run_module(job_count=num_jobs)
    except Exception as e:
        shell.clear_working_directory()
        raise e
    shell.check_module_output()
    shell.write_sample_sheet()
    shell.clear_working_directory()
    if shell.pack_output: shell.fold_output()
