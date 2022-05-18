from ardetype_utilities import parse_folder, create_sample_sheet, read_config, edit_config, write_config, submit_module_job, run_module_cluster, check_job_completion, check_module_output, edit_sample_sheet, validate_config
import os, warnings, re
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

class Module():
    '''Class represents single module of the ardetype pipeline, that can interact with other modules through'''
    modules = {
        "core" : {
            "targets":[
                "_contigs.fasta",
                "_bact_reads_classified_1.fastq.gz", 
                "_bact_reads_classified_2.fastq.gz",
                "_bact_reads_unclassified_1.fastq.gz",
                "_bact_reads_unclassified_2.fastq.gz",
                "_kraken2_contigs_report.txt",
                "_kraken2_host_filtering_report.txt"
            ],
            "patterns":{
                "inputs":[".fastq.gz"],
                "sample_sheet":"(_R[1,2]_001.fastq.gz|_[1,2].fastq.gz)"
            },
            "job_name":"bact_core"
        },
        "shell" : {
            "targets":[
                ".rgi.txt",
                ".rgi.json",
                "_mlst_output.csv",
                "_resfinder/pheno_table.txt",
                "_resfinder/ResFinder_Hit_in_genome_seq.fsa",
                "_resfinder/ResFinder_Resistance_gene_seq.fsa",
                "_resfinder/ResFinder_results_tab.txt",
                "_resfinder/ResFinder_results.txt",
                "_amrpp/ResistomeResults/AMR_analytic_matrix.csv"
            ],
            "patterns":{
                "inputs":[".fastq.gz", "_contigs.fasta"],
                "sample_sheet":"(_R[1,2]_001.fastq.gz|_[1,2].fastq.gz|_contigs.fasta)"
            },
            "job_name":"bact_shell"
        }
    }

    def __init__(self, module_name, input_path, module_config, output_path, run_mode) -> None:
        self.run_mode = run_mode
        self.job_id = None
        self.module_name = module_name
        self.input_path = input_path
        self.output_path = output_path
        self.target_list = None
        self.sample_sheet = None
        self.config_file_path = f'{os.path.abspath(self.output_path)}/config.yaml'
        self.cluster_config_path = 'cluster.yaml'
        self.config_file = read_config(module_config) if isinstance(module_config, str) else module_config
        self.input_dict = {}
        self.sample_sheet = None

    def fill_input_dict(self):
        '''Methods fills self.input_dict using self.input_path and self.module_name by
        mapping each file format to the list of files of that format, found in the self.input_path'''
        for format in Module.modules[self.module_name]['patterns']['inputs']:
            self.input_dict[format] = parse_folder(self.input_path, format)
   
    def fill_sample_sheet(self):
        '''
        Method initializes self.sample_sheet to pandas dataframe, using self.input_dict and self.module_name (restricted to fastq & fasta inputs)
        '''
        if len(self.input_dict) < 2:
            self.sample_sheet = create_sample_sheet(self.input_dict[".fastq.gz"],Module.modules[self.module_name]['patterns']['sample_sheet'],mode=0)
        else:
            self.sample_sheet = create_sample_sheet(self.input_dict[".fastq.gz"],Module.modules[self.module_name]['patterns']['sample_sheet'],mode=0)
            fasta_dict = {re.sub("_contigs.fasta","",os.path.basename(contig)):contig for contig in self.input_dict["_contigs.fasta"]}
            self.sample_sheet = edit_sample_sheet(self.sample_sheet,fasta_dict,'fa')

    def add_fasta_samples(self):
        '''Adds fa column with fasta files to the self.sample_sheet dataframe'''
        fasta_dict = {re.sub("_contigs.fasta","",os.path.basename(contig)):contig for contig in self.input_dict["_contigs.fasta"]}
        self.sample_sheet = edit_sample_sheet(self.sample_sheet,fasta_dict,'fa')


    def fill_target_list(self):
        '''Method fills self.target_list using data stored in self.sample_sheet instance variable'''
        self.target_list = [f'{self.output_path}{id}{tmpl}' for id in self.sample_sheet['sample_id'] for tmpl in Module.modules[self.module_name]['targets']]


    def make_output_dir(self):
        '''Creates output directory (if not present in the file system) using self.output_path'''
        if not os.path.exists(self.output_path): os.makedirs(self.output_path)


    def write_sample_sheet(self):
        '''Creates sample_sheet.csv file in the self.output_path folder, using self.sample_sheet'''
        self.sample_sheet.to_csv(f"{self.output_path}sample_sheet.csv", header=True, index=False)


    def add_module_targets(self):
        '''Updates self.config_file, using self.module_name'''
        output_code = edit_config(self.config_file, f"{self.module_name}_target_files", self.target_list)
        validation_code = validate_config(self.config_file)
        if not output_code == 0: raise Exception(f'Config editing failed with error code {output_code}')
        elif not validation_code == 0: raise Exception(f'Config validation failed with error code {validation_code}')


    def add_output_dir(self):
        '''Updates self.config_file using self.output_path'''
        output_code = edit_config(self.config_file, "output_directory", self.output_path)
        validation_code = validate_config(self.config_file)
        if not output_code == 0: raise Exception(f'Config editing failed with error code {output_code}')
        elif not validation_code == 0: raise Exception(f'Config validation failed with error code {validation_code}')


    def write_module_config(self):
        '''Writes self.config_file to the self.output_path'''
        write_config(self.config_file, f'{self.output_path}config.yaml')


    def run_module(self):
        '''Runs module on hpc as job or as snakemake submitter (on login node), based on self.run_mode value (True - job, False - submitter)'''
        if self.run_mode:
            self.job_id = submit_module_job(self.module_name,self.config_file_path, self.output_path)
            check_job_completion(self.job_id,Module.modules[self.module_name]['job_name'],sleeping_time=5,output_dir=self.output_path)
        else:
            run_module_cluster(self.module_name,self.config_file_path, self.cluster_config_path, 12)


    def check_module_output(self):
        '''Checks if output files are generated according to self.module_name and adds check_note_{self.module_name} column 
        to the self.sample_sheet dataframe, where boolean value is stored for each expected file'''
        check_dict = check_module_output(file_list=self.target_list)
        id_check_dict = {id:"" for id in self.sample_sheet['sample_id']}
        for file in check_dict:
            two_dirs_up = os.path.basename(os.path.dirname(os.path.dirname(file)))+"/"+os.path.basename(os.path.dirname(file))+"/"+os.path.basename(file)
            id = os.path.basename(re.sub("("+"|".join(Module.modules[self.module_name]['targets'])+")","",two_dirs_up))
            id_check_dict[id] += f"|{file}:{check_dict[file]}"
        self.sample_sheet = edit_sample_sheet(self.sample_sheet, id_check_dict, f"check_note_{self.module_name}")


    def supply_sample_sheet(self):
        '''Returns self.sample_sheet dataframe object'''
        return self.sample_sheet


    def receive_sample_sheet(self, sample_sheet):
        '''Inializes self.sample_sheet with external sample_sheet dataframe (used to connect modules)'''
        self.sample_sheet = sample_sheet
