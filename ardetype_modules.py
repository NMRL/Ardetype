from ardetype_utilities import *
import os, warnings, re, sys, subprocess, shutil, time
from itertools import chain
from getpass import getuser

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

class Module:
    '''Class represents single module of the ardetype pipeline'''

    def __init__(self, module_name: str, input_path: str, module_config, output_path: str, run_mode: bool, job_name: str, patterns: dict, targets: list, requests: dict, snakefile_path: str, cluster_config_path: str) -> None:
        self.run_mode = run_mode
        self.job_id = None
        self.taxonomy_dict = None
        self.module_name = module_name
        self.input_path = input_path
        self.output_path = output_path
        self.target_list = None
        self.sample_sheet = None
        self.aggr_taxonomy_path = f'{os.path.abspath(self.output_path)}/{self.module_name}_aggregated_taxonomy.json'
        self.config_file_path = f'{os.path.abspath(self.output_path)}/config.yaml'
        self.cluster_config_path = cluster_config_path
        self.config_file = read_yaml(module_config) if isinstance(module_config, str) else module_config
        self.input_dict = {}
        self.sample_sheet = None
        self.patterns = patterns
        self.job_name = job_name
        self.targets = targets
        self.requests = requests
        self.snakefile_path = snakefile_path


    def fill_input_dict(self, substring_list=['reads_unclassified', 'reads_classified']):
        '''Fills self.input_dict using self.input_path and self.module_name by
        mapping each file format to the list of files of that format, found in the self.input_path, excluding files that contain substrings in their names (supply None to avoid excluding files)'''
        for format in self.patterns['inputs']:
            self.input_dict[format] = parse_folder(self.input_path,substr_lst=substring_list, file_fmt_str=format)
   

    def fill_sample_sheet(self):
        '''
        Method initializes self.sample_sheet to pandas dataframe, using self.input_dict and self.module_name (restricted to fastq & fasta inputs)
        '''
        if len(self.input_dict) < 2:
            self.sample_sheet = create_sample_sheet(self.input_dict[".fastq.gz"],self.patterns['sample_sheet'],mode=0)
        else:
            self.sample_sheet = create_sample_sheet(self.input_dict[".fastq.gz"],self.patterns['sample_sheet'],mode=0)
            fasta_dict = {re.sub("_contigs.fasta","",os.path.basename(contig)):contig for contig in self.input_dict["_contigs.fasta"]}
            self.sample_sheet = map_new_column(self.sample_sheet,fasta_dict,'sample_id','fa')


    def fill_target_list(self, taxonomy_based=False):
        '''Method fills self.target_list using data stored in self.sample_sheet instance variable'''
        if taxonomy_based:
            self.target_list = [f'{self.output_path}{id}{tmpl}' for idx, id in enumerate(self.sample_sheet['sample_id']) for tmpl in self.targets[self.sample_sheet['taxonomy'][idx]]]
        else:
            self.target_list = [f'{self.output_path}{id}{tmpl}' for id in self.sample_sheet['sample_id'] for tmpl in self.targets]
            


    def make_output_dir(self):
        '''Creates output directory (if not present in the file system) using self.output_path'''
        if not os.path.exists(self.output_path): os.makedirs(self.output_path)


    def write_sample_sheet(self):
        '''Creates sample_sheet.csv file in the self.output_path folder, using self.sample_sheet'''
        self.sample_sheet.to_csv(f"{self.output_path}sample_sheet.csv", header=True, index=False)


    def add_module_targets(self):
        '''Updates self.config_file, using self.module_name'''
        output_code = edit_nested_dict(self.config_file, f"{self.module_name}_target_files", self.target_list)
        validation_code = validate_yaml(self.config_file)
        if not output_code == 0: raise Exception(f'Config editing failed with error code {output_code}')
        elif not validation_code == 0: raise Exception(f'Config validation failed with error code {validation_code}')


    def add_output_dir(self):
        '''Updates self.config_file using self.output_path'''
        output_code = edit_nested_dict(self.config_file, "output_directory", self.output_path)
        validation_code = validate_yaml(self.config_file)
        if not output_code == 0: raise Exception(f'Config editing failed with error code {output_code}')
        elif not validation_code == 0: raise Exception(f'Config validation failed with error code {validation_code}')


    def write_module_config(self):
        '''Writes self.config_file to the self.output_path'''
        write_yaml(self.config_file, f'{self.output_path}config.yaml')


    def check_module_output(self):
        '''Checks if output files are generated according to self.module_name and adds check_note_{self.module_name} column 
        to the self.sample_sheet dataframe, where boolean value is stored for each expected file'''
        check_dict = check_file_existance(file_list=self.target_list)
        id_check_dict = {id:"" for id in self.sample_sheet['sample_id']}
        for file in check_dict:
            two_dirs_up = os.path.basename(os.path.dirname(os.path.dirname(file)))+"/"+os.path.basename(os.path.dirname(file))+"/"+os.path.basename(file)
            if isinstance(self.targets, list):
                id = os.path.basename(re.sub("("+"|".join(self.targets)+")","",two_dirs_up))
            elif isinstance(self.targets, dict):
                id = os.path.basename(re.sub("("+"|".join(chain.from_iterable(self.targets.values()))+")","",two_dirs_up))
            id_check_dict[id] += f"|{file}:{check_dict[file]}"
        self.sample_sheet = map_new_column(self.sample_sheet, id_check_dict, 'sample_id', f"check_note_{self.module_name}")


    def supply_sample_sheet(self):
        '''Returns self.sample_sheet dataframe object'''
        return self.sample_sheet


    def receive_sample_sheet(self, sample_sheet):
        '''Inializes self.sample_sheet with external sample_sheet dataframe (used to connect modules)'''
        self.sample_sheet = sample_sheet


    def add_fasta_samples(self):
        '''Adds fa column with _contigs.fasta files to the self.sample_sheet dataframe'''
        fasta_dict = {re.sub("_contigs.fasta","",os.path.basename(contig)):contig for contig in self.input_dict["_contigs.fasta"]}
        self.sample_sheet = map_new_column(self.sample_sheet,fasta_dict,'sample_id','fa')


    def remove_invalid_samples(self, connect_from_module_name):
        '''
        Given supplier module name, removes samples that lack files, required by the current module.
        If all samples are removed, returns 1 (int).
        '''
        if self.requests['check'] is not None:
            if isinstance(self.requests['check'],str): 
                self.sample_sheet = self.sample_sheet[self.sample_sheet[f'check_note_{connect_from_module_name}'].str.contains(self.requests['check'])].reset_index(drop=True)
            else:
                for request in self.requests['check']:
                    self.sample_sheet = self.sample_sheet[self.sample_sheet[f'check_note_{connect_from_module_name}'].str.contains(request)].reset_index(drop=True)
        if self.requests['taxonomy'] is not None:
            self.sample_sheet = self.sample_sheet[self.sample_sheet['taxonomy'].str.contains("("+"|".join(self.requests['taxonomy'])+")")].reset_index(drop=True)
        if self.sample_sheet.empty:
            return 1


    def submit_module_job(self, jobscript_path):
        """
        Given path to the job_script perform job submition to RTU HPC cluster, setting self.job_id to the bytestring representing job_id.
        """
        shutil.copy(jobscript_path, f'{self.output_path}ardetype_jobscript.sh')

        try:
            self.job_id = subprocess.check_output(['qsub', '-F', f'{self.snakefile_path} {self.config_file_path}', f'{self.output_path}ardetype_jobscript.sh'])
            os.system(f"rm {self.output_path}ardetype_jobscript.sh")
        except subprocess.CalledProcessError as msg:
            raise Exception(f"Job submission error: {msg}")
        

    def check_job_completion(self, sleeping_time=5):
        """
        Given sleeping time (in seconds) between checks (int), checks if the job is complete every n seconds. 
        Waiting is finished when the job status is C (complete), when the file with the job standard output is moved to self.output_path.
        """
        self.job_id = self.job_id.decode('UTF-8').strip()
        search_string = f"{self.job_id}.*{getuser()}.*[RCQ]"
        print(f'Going to sleep until ardetype/{self.module_name}/{self.job_id} job is finished')
        while True:
            qstat = os.popen("qstat").read()
            check_job = re.search(search_string,qstat).group(0)
            print(f"{check_job} : {time.ctime(time.time())}")
            if check_job[-1] == "C":
                print(f'Finished waiting: ardetype/{self.module_name}/{self.job_id} is complete')
                break
            time.sleep(sleeping_time)
            job_report = f"*o{self.job_id.split('.')[0]}"
            os.system(f"mv {job_report} {self.output_path}/{self.module_name}_{self.job_name}_{self.job_id}.txt")


    def run_module_cluster(self, job_count):
        '''
        Given number of jobs to run in parallel (int) runs module on the login node of the HPC cluster, 
        allowing the snakemake to do job submissions to the computing nodes automatically.     
        '''
 
        qsub_command = '"qsub -N {cluster.jobname} -l nodes={cluster.nodes}:ppn={cluster.ppn},pmem={cluster.pmem},walltime={cluster.walltime} -q {cluster.queue} -j {cluster.jobout} -o {cluster.outdir} -V"'
        shell_command = f'''
        eval "$(conda shell.bash hook)";
        conda activate /mnt/home/$(whoami)/.conda/envs/mamba_env/envs/snakemake; 
        snakemake --jobs {job_count} --cluster-config {self.cluster_config_path} --cluster-cancel qdel --configfile {self.config_file_path} --snakefile {self.snakefile_path} --keep-going --use-envmodules --use-conda --conda-frontend conda --rerun-incomplete --latency-wait 30 --cluster {qsub_command}'''
        try:
            subprocess.check_call(shell_command, shell=True)
        except subprocess.CalledProcessError as msg:
            sys.exit(f"Module process running error: {msg}")


    def run_module(self, job_count, jobscript_path='./ardetype_jobscript.sh'):
        '''Runs module on hpc as job or as snakemake submitter (on login node), based on self.run_mode value (True - job, False - submitter)'''
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
