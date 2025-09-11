from subscripts.ardetype_utilities import Ardetype_housekeeper as hk
import os, warnings, sys, subprocess as sp, datetime, pandas as pd, re
from pathlib import Path
from shutil import move
sys.path.insert(0, os.path.dirname(Path(__file__).absolute()))

from src.modules import Module

#Suppressing pandas warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)


#Reading data used to build module objects
ardetype_path = os.path.dirname(Path(__file__).parents[0].absolute())
module_data   = hk.read_json_dict(os.path.join(ardetype_path,'config_files/json/module_data.json'))

###############################################################
# Wrapper to extract tool versions
###############################################################
class Wrapper():
    '''Toolkit class to store wrapper methods & attributes for different tools'''

    #pipeline configuration saved at module import
    def __init__(self):
        ardetype_path = str(Path(os.path.abspath('./')))
        self._config_dict   = hk.read_yaml(os.path.join(ardetype_path,'config_files/yaml/config_modular_local.yaml'))
        self._tool_ref_map  = hk.read_json_dict(os.path.join(ardetype_path, 'config_files/json/specific_tool_map.json'))
        self._db_vers_map = None
        self._tool_vers_map = None

        #Updating paths in config to match user file system
        for p in ["databases"]:
            if not os.path.isabs(self._config_dict[p]):
                self._config_dict[p] = os.path.join(ardetype_path, self._config_dict[p])

        for substr in ["_database", "_db", "_sif"]:
            hk.join_sif_paths(self._config_dict, substr, ardetype_path)

    def _get_datestamp(self, common_dir:str, name:str)->str:
        '''
        Returns date of latest file/folder manipulation
        given path to common parent folder and specific file name.
        '''
        datestamp = datetime.datetime.fromtimestamp(os.path.getmtime(os.path.join(common_dir, name))).strftime('%Y-%m-%d')
        return datestamp


    def set_db_vers_map(self):
        '''Setter for db_vers_map'''
        # Create base path for databases from the configuration
        db_path = self._config_dict["databases"]

        # Refactor the original dictionary with variables and dictionary comprehensions

        self._db_vers_map = {
            "kraken2": {
                "filter_host": self._get_datestamp(f'{db_path}db-kraken2/', os.path.basename(os.path.dirname(self._config_dict['core_tool_configs']['kraken2']['human_db']))),
                "classify_reads": self._get_datestamp(f'{db_path}db-kraken2/', os.path.basename(os.path.dirname(self._config_dict['core_tool_configs']['kraken2']['bact_db']))),
                "classify_contigs": self._get_datestamp(f'{db_path}db-kraken2/', os.path.basename(os.path.dirname(self._config_dict['core_tool_configs']['kraken2']['bact_db']))),
            },
            "chewbbaca": {
                k: self._get_datestamp(f'{db_path}{self._config_dict["tip_tool_configs"]["chewbbaca"]["chewbbaca_dtb"]}', v)
                for k, v in {
                    "Salmonella enterica" : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['salmonella']['db'],
                    "Escherichia coli" : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['escherichia']['db'],
                    "Enterobacter"  : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['escherichia']['db'],
                    "Shigella"  : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['escherichia']['db'],
                    "Streptococcus pyogenes" : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['streptococcus pyogenes']['db'],
                    "Streptococcus pneumoniae" : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['streptococcus pneumoniae']['db'],
                    "Acinetobacter baumannii" : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['baumannii']['db'],
                    "Campylobacter"  : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['campylobacter']['db'],
                    "Klebsiella"  : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['klebsiella']['db'],
                    "Legionella pneumophila" : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['legionella']['db'],
                    "Listeria monocytogenes" : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['listeria']['db'],
                    "Staphylococcus aureus" : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['staphylococcus aureus']['db'],
                    "Enterococcus faecalis" : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['enterococcus faecalis']['db'],
                    "Enterococcus faecium" : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['enterococcus faecium']['db'],
                    "Pseudomonas aeruginosa" : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['pseudomona']['db'],
                    "Yersinia"  : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['yersinia']['db'], 
                    "Clostridioides"  : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['clostridioides']['db'],
                    "Clostridium"  : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['clostridium']['db'],
                    "Mycobacterium"  : self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']['mycobacterium']['db'],
                }.items()
            },
            "resfinder": self._get_datestamp(self._config_dict['shell_tool_configs']['resfinder']['resfinder_db'], 'resfinder_db'),
            "plasmidfinder": self._get_datestamp(f'{db_path}', 'plasmidfinder_db_new'),
            "virulencefinder": self._get_datestamp(f'{db_path}', 'virulencefinder_db_new'),
            "plasmidtype_abaumannii": self._get_datestamp(f'{db_path}', 'acinetobacterplasmidtyper.fasta')
        }

    def set_tool_vers_map(self):
        '''Setter for tool_vers_map'''
        self._tool_vers_map = {
            "plasmidfinder": 
                datetime.datetime.fromtimestamp(os.path.getmtime(self._config_dict["shell_tool_configs"]["plasmidfinder"]["plasmidfinder_sif"])).strftime('%Y-%m-%d'),
            "resfinder":
                sp.run(
                f'singularity --silent run {self._config_dict["resfinder_sif"]} python -m resfinder --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip(),
            "virulencefinder":
                sp.run(f'singularity --silent run {self._config_dict["shell_tool_configs"]["virulencefinder"]["virulencefinder_sif"]} python -m virulencefinder --version 2> /dev/null', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip(),
            "quast":
                sp.run(
                f'singularity --silent run {self._config_dict["quast_sif"]} quast --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split('\n')[-1],
            "rgi":
                sp.run(
                f'singularity --silent run {self._config_dict["rgi_sif"]} rgi main --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip(),
            "kraken2":
                sp.run(
                f'singularity --silent run {self._config_dict["kraken2_sif"]} kraken2 --version',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split('\n')[0],
            "amrfinder+":
                sp.run(f'singularity --silent run {self._config_dict["amrfinderplus_sif"]} amrfinder --version',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip(),
            "fastp": 
                sp.run(f'singularity --silent run {self._config_dict["fastp_sif"]} fastp --version',
                stderr=sp.PIPE, shell=True).stderr.decode('utf-8').strip(),
            "mob-suite":
                sp.run(f'singularity --silent run {self._config_dict["mob_suite_sif"]} mob_typer --version',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').split(" ")[-1].strip(),
            "mlst":
                sp.run(f'singularity --silent run {self._config_dict["mlst_sif"]} mlst --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "shovill":
                sp.run(f'singularity --silent run {self._config_dict["shovill_sif"]} shovill --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "meningotype":
                sp.run(f'singularity --silent run {self._config_dict["meningotype_nmeningitidis_sif"]} meningotype --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "legsta":
                sp.run(f'singularity --silent run {self._config_dict["legsta_lpneumophila_sif"]} legsta --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "hicap":
                sp.run(f'singularity --silent run {self._config_dict["hicap_hinfluenzae_sif"]} hicap --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "kleborate":
                sp.run(f'singularity --silent run {self._config_dict["kleborate_kpneumoniae_sif"]} kleborate --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "agrvate":
                sp.run(f'singularity --silent run {self._config_dict["agrvate_saureus_sif"]} agrvate --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "spatyper":
                sp.run(f'singularity --silent run {self._config_dict["spatyper_saureus_sif"]} spaTyper --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "staphopia-sccmec":
                sp.run(f'singularity --silent run {self._config_dict["sccmec_saureus_sif"]} staphopia-sccmec --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "emmtyper":
                sp.run(f'singularity --silent run {self._config_dict["emmtyper_spyogenes_sif"]} emmtyper --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "lissero":
                sp.run(f'singularity --silent run {self._config_dict["lissero_lmonocytogenes_sif"]} lissero --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "sistr":
                sp.run(f'singularity --silent run {self._config_dict["sistr_senterica_sif"]} sistr --version ',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "seqsero2":
                sp.run(f'singularity --silent run {self._config_dict["seqsero2_senterica_sif"]} SeqSero2_package.py --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "ectyper":
                sp.run(f'singularity --silent run {self._config_dict["ectyper_ecoli_sif"]} ectyper --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip(),
            "stecfinder":
                sp.run(f'singularity --silent run {self._config_dict["stecfinder_ecoli_sif"]} stecfinder --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "seroba":
                sp.run(f'singularity --silent run {self._config_dict["seroba_spneumoniae_sif"]} seroba version',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip(),
            "chewbbaca":
                sp.run(f'singularity --silent run {self._config_dict["tip_tool_configs"]["chewbbaca"]["chewbbaca_sif"]} chewBBACA.py --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            'lrefinder':
                sp.run(f'singularity --silent run  {self._config_dict["lrefinder_efaecium_efaecalis_sif"]} LRE-Finder.py -v', 
                stderr=sp.PIPE, shell=True).stderr.decode('utf-8').strip().split('\n')[-1],
            "shigatyper":
                sp.run(f'singularity --silent run {self._config_dict["shigatyper_shigella_sif"]} shigatyper --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "flye":
                sp.run(f'singularity --silent run {self._config_dict["flye_sif"]} flye --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "snikt":
                sp.run(f'singularity --silent run {self._config_dict["snikt_sif"]} snikt.R --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "filtlong":
                sp.run(f'singularity --silent run {self._config_dict["filtlong_sif"]} filtlong --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "circlator":
                sp.run(f'singularity --silent run {self._config_dict["circlator_sif"]} circlator version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "polypolish":
                sp.run(f'singularity --silent run {self._config_dict["polypolish_sif"]} polypolish --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "polca":"4.1.0",
            "medaka": sp.run(f'singularity --silent run {self._config_dict["medaka_sif"]} medaka --version 2> /dev/null', stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split('\n')[-1].split()[1],
            "capybara":"1.0",
            "plasmidtype_abaumannii":"2022",
            "hybracter": "0.11.0"
        }


    def report_tool_versions(self, output_path:str, sample_ids:list) -> None:
        '''Aggregates tool versions, references and database versions in single pandas dataframe and generates a csv file named software_log.csv in output_path'''
        specific_ref = {self._tool_ref_map[k]['tool']:self._tool_ref_map[k]['link'] for k in self._tool_ref_map if k != 'agnostic'}
        agnostic_ref = {self._tool_ref_map['agnostic'][k]['tool'] : self._tool_ref_map['agnostic'][k]['link'] for k in self._tool_ref_map['agnostic']}
        reference_mapa = {self._tool_ref_map['agnostic'][k]['tool'] : self._tool_ref_map['agnostic'][k]['reference'] for k in self._tool_ref_map['agnostic']}
        reference_maps = {self._tool_ref_map[k]['tool']:self._tool_ref_map[k]['reference'] for k in self._tool_ref_map if k != 'agnostic'}
        total_ref = {**specific_ref, **agnostic_ref}
        total_cit = {**reference_mapa, **reference_maps}

        df = pd.DataFrame.from_dict({'tool':[t for t in self._tool_vers_map], 'tool_version':[v for v in self._tool_vers_map.values()]})
        df['link'] = df['tool'].map(total_ref)
        df['reference'] = df['tool'].map(total_cit)
        df['run_parameters'] = ['default' for _ in df.index]

        db_map = {}
        for t in self._db_vers_map:
            if t != 'chewbbaca':
                if isinstance(self._db_vers_map[t], dict):
                    db_map[t] = " | ".join([f"{r} : {v}" for r,v in self._db_vers_map[t].items()])             
                else:
                    db_map[t] = self._db_vers_map[t]
        df['db_versions'] = df['tool'].map(db_map)
        df['db_versions'] = df['db_versions'].fillna('not_applicable')
        df.insert(0, 'analysis_batch_id', [os.path.basename(os.path.dirname(output_path)) for _ in df.index])
        df['ardetype_version'] = self._config_dict['ardetype_version']

        #Add cgmlst databases as tools
        for org in self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs']:
            try:
                timestamp = [self._db_vers_map['chewbbaca'][o] for o in self._db_vers_map['chewbbaca'] if org.lower().replace(' ', '_') in o.lower().replace(' ', '_')][0]
            except IndexError:
                timestamp = 'NA'
            df.loc[len(df)] = {
                "analysis_batch_id":os.path.basename(os.path.dirname(output_path)),
                "tool":f"cgmlst_{org.lower()}",
                "tool_version":"not_applicable",
                "run_parameters":"default",
                "db_versions":timestamp,
                "reference":self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs'][org]['reference'],
                "link":self._config_dict['tip_tool_configs']['chewbbaca']['chewbbaca_orgs'][org]['origin'],
                "ardetype_version":self._config_dict['ardetype_version'],
                }

        df = df[["analysis_batch_id","tool","tool_version","run_parameters","db_versions","ardetype_version","reference","link"]]
        log_path = os.path.join(output_path, 'software_log.csv')
        # if not os.path.isfile(log_path):
        df.to_csv(log_path, header=True, index=False)
        for sid in sample_ids:
            df.to_csv(os.path.join(output_path, f'{sid}_tool_log.csv'), header=True, index=False)


###############################################################
# Extending base Module class to fit the needs of the pipeline
###############################################################

class Ardetype_module(Module):
    '''Class extends Module and implements pipeline-specific methods'''
    def __init__(self, *args, **kwargs):
        '''
        Method is overridden because post-argument parsing validation of
        command-line arguments is required by some modules but not others, 
        which cannot be covered by altering the configuration file.
        '''
        if kwargs['input_path'] is None:
            sys.exit(f'Input must be included to run the pipeline in any mode except log_analysis and merge.')
        super(Ardetype_module, self).__init__(*args, **kwargs) #running method as it is defined in the base class
        self.job_log_path = os.path.join(ardetype_path,'ardetype_job_logs/')
        self.status_script = self.config_file['status_script_path']
        self.unfold_output()

        #if reprocess
        if self.force_all:
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            #change timestamp
            if re.match(r'_[0-9]{8}_[0-9]{6}/', self.output_path[-17:]):
                new_path = self.output_path[:-17] + f"_{timestamp}/"
            else:
                #if timestamp is missing - add timestamp
                new_path = os.path.dirname(self.output_path)+f"_{timestamp}/"
            try:
                move(self.output_path, new_path)
            except FileNotFoundError as e:
                print(f'Was unable to find {self.output_path} directory - please check if the timestamp is correct.', file=sys.stderr)
                sys.exit(1)
            if new_path not in self.aggr_taxonomy_path:                
                self.aggr_taxonomy_path  = os.path.join(new_path, os.path.basename(self.aggr_taxonomy_path))
            if new_path not in self.config_file_path:
                self.config_file_path    = os.path.join(new_path, os.path.basename(self.config_file_path))
            
            if os.path.abspath(self.output_path) == os.path.abspath(self.input_path):
                self.input_path = os.path.abspath(new_path)

            self.output_path = new_path
            

        #if not reprocess but no timestamp
        elif not re.match(r'_[0-9]{8}_[0-9]{6}/', self.output_path[-17:]):
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            new_path  = os.path.dirname(self.output_path)+f"_{timestamp}/"
            if os.path.isdir(self.output_path):
                move(self.output_path, new_path)
            if new_path not in self.aggr_taxonomy_path:
                self.aggr_taxonomy_path  = os.path.join(new_path, os.path.basename(self.aggr_taxonomy_path))
            if new_path not in self.config_file_path:
                self.config_file_path    = os.path.join(new_path, os.path.basename(self.config_file_path))
            if os.path.abspath(self.output_path) == os.path.abspath(self.input_path):
                self.input_path = os.path.abspath(new_path)
            self.output_path    = new_path
           

            
            


    def config_cluster(self) -> None:
        '''
        Change pre-set parameters in default config_cluster.yaml,
        copy cluster.yaml into self.output_path/logs/
        and change self.config_file_path to self.output_path/logs/
        '''
        self.job_log_path     = os.path.join(self.output_path,'logs/')
        cluster_config        = hk.read_yaml(self.cluster_config_path)
        os.makedirs(self.job_log_path , exist_ok=True)
        hk.edit_nested_dict(cluster_config,'outdir', self.job_log_path)
        if not self.run_local:
            hk.write_yaml(cluster_config, os.path.join(self.output_path,'config_cluster.yaml'))
        self.cluster_config_path = os.path.join(self.output_path,'config_cluster.yaml')


###############################################
# Defining wrapper functions to call from main
###############################################

def run_all(args, num_jobs):
    '''Wrapper function to run all modules sequentially.'''
    core = Ardetype_module(
            module_name='core',
            input_path          = args.input,
            module_config       = args.config,
            output_path         = args.output_dir,
            run_mode            = args.submit_modules,
            run_local           = not args.hpc,
            force_all           = args.force_all,
            snakemake_cpus      = args.max_cpus,
            pack_output         = not args.skip_packing,
            retry_times         = args.retry_times,
            rules_to_rerun      = args.force_rules,
            job_name            = module_data['core']['job_name'],
            patterns            = module_data['core']['patterns'],
            targets             = module_data['core']['targets'],
            requests            = module_data['core']['requests'],
            snakefile_path      = module_data['snakefiles']['core'],
            cluster_config_path = module_data['cluster_config'],
            nanopore_mode       = args.nanopore_only,
            fasta_mode          = args.fasta
            )

    shell = Ardetype_module(
        module_name         = 'shell', 
        input_path          = core.output_path if core.output_path == core.input_path else core.input_path, 
        module_config       = core.config_file, 
        output_path         = core.output_path, 
        run_mode            = args.submit_modules,
        run_local           = not args.hpc,
        force_all           = args.force_all,
        snakemake_cpus      = args.max_cpus,
        pack_output         = not args.skip_packing,
        retry_times         = args.retry_times,
        rules_to_rerun      = args.force_rules,
        job_name            = module_data['shell']['job_name'],
        patterns            = module_data['shell']['patterns'],
        targets             = module_data['shell']['targets'],
        requests            = module_data['shell']['requests'],
        snakefile_path      = module_data['snakefiles']['shell'],
        cluster_config_path = module_data['cluster_config'],
        nanopore_mode       = args.nanopore_only,
        fasta_mode          = args.fasta
        )
    tip = Ardetype_module(
        module_name         = 'tip', 
        input_path          = core.output_path if core.output_path == core.input_path else core.input_path, 
        module_config       = shell.config_file, 
        output_path         = core.output_path, 
        run_mode            = args.submit_modules,
        run_local           = not args.hpc,
        force_all           = args.force_all,
        snakemake_cpus      = args.max_cpus,
        pack_output         = not args.skip_packing,
        retry_times         = args.retry_times,
        rules_to_rerun      = args.force_rules,
        job_name            = module_data['tip']['job_name'],
        patterns            = module_data['tip']['patterns'],
        targets             = module_data['tip']['targets'],
        requests            = module_data['tip']['requests'],
        snakefile_path      = module_data['snakefiles']['tip'],
        cluster_config_path = module_data['cluster_config'],
        nanopore_mode       = args.nanopore_only,
        fasta_mode          = args.fasta
    )
    shape = Ardetype_module(
        module_name         = 'shape',
        input_path          = core.output_path if core.output_path == core.input_path else core.input_path, 
        module_config       = tip.config_file,
        output_path         = core.output_path,
        run_mode            = args.submit_modules,
        run_local           = not args.hpc,
        force_all           = args.force_all,
        snakemake_cpus      = args.max_cpus,
        pack_output         = not args.skip_packing,
        retry_times         = args.retry_times,
        rules_to_rerun      = args.force_rules,
        job_name            = module_data['shape']['job_name'],
        patterns            = module_data['shape']['patterns'],
        targets             = module_data['shape']['targets'],
        requests            = module_data['shape']['requests'],
        snakefile_path      = module_data['snakefiles']['shape'],
        cluster_config_path = module_data['cluster_config'],
        nanopore_mode       = args.nanopore_only,
        fasta_mode          = args.fasta
    )

    #Running core
    print(f'{core.module_name}: Checking input', file=sys.stdout)
    if args.nanopore_only:
        print(f'{core.module_name}: Nanopore mode requested - configuring', file=sys.stdout)
        core.snakefile_path = module_data['snakefiles']['core_ont']
        core.fill_input_dict(pattern_path='ONT')
        core.fill_sample_sheet(pattern=core.patterns['inputs']['ONT'])
    elif args.fasta:
        print(f'{core.module_name}: Fasta mode requested - configuring', file=sys.stdout)
        core.snakefile_path = module_data['snakefiles']['core_fasta']
        core.fill_input_dict(pattern_path='FASTA')
        core.fill_sample_sheet(pattern=core.patterns['inputs']['FASTA'])
    else:
        core.fill_input_dict(pattern_path='ILL')
        core.fill_sample_sheet(pattern=core.patterns['inputs']['ILL'])
    print(f'{core.module_name}: Adding targets to config file', file=sys.stdout)
    core.make_output_dir()
    core.get_sample_groups(fasta=args.fasta)
    core.write_sample_sheet()
    core.fill_target_list(grouped=True)
    core.add_module_targets()
    core.add_output_dir()
    core.config_cluster()
    core.write_module_config()
    print(f'{core.module_name}: Launching snakemake', file=sys.stdout)
    core.run_module(job_count=num_jobs)
    print(f'{core.module_name}: Verifying output completion', file=sys.stdout)
    core.check_module_output()
    core.add_taxonomy_column()
    core.write_sample_sheet()
    print(f'{core.module_name}: Finished', file=sys.stdout)

  
    #Connecting core to shell
    print(f'{shell.module_name}: Parsing core output', file=sys.stdout)
    shell.receive_sample_sheet(core.supply_sample_sheet())
    print(f'{shell.module_name}: Removing samples lacking required files', file=sys.stdout)
    samples_cleared = shell.remove_invalid_samples(connect_from_module_name='core') #in dry run mode none of the rules are executed, hence all samples will be removed, causing error
    shell.save_removed()
    if samples_cleared == 1: 
        if core.pack_output: core.fold_output()
        print(f'{shell.module_name}: Missing files requested by bact_shell', file=sys.stderr)
        sys.exit(1)

    #Running shell
    if args.nanopore_only or args.fasta:
        print(f'{shell.module_name}: {"Nanopore" if args.nanopore_only else "Fasta"} mode requested - configuring', file=sys.stdout)
        shell.snakefile_path = module_data['snakefiles']['shell_ont']
        shell.fill_input_dict(pattern_path='ONT')
    else:
        shell.fill_input_dict(pattern_path='ILL/FUL')

    print(f'{shell.module_name}: Adding targets to config file', file=sys.stdout)
    shell.add_fasta_samples()
    shell.write_sample_sheet()
    shell.fill_target_list()
    shell.add_module_targets()
    shell.config_cluster()
    shell.write_module_config()

    print(f'{shell.module_name}: Launching snakemake', file=sys.stdout)
    shell.run_module(job_count=num_jobs)

    print(f'{shell.module_name}: Verifying output completion', file=sys.stdout)
    shell.check_module_output()
    shell.add_taxonomy_column()
    shell.write_sample_sheet()
    print(f'{shell.module_name}: Finished', file=sys.stdout)

    # Connecting shell & core to tip/shape
    print(f'{tip.module_name}: Checking required input', file=sys.stdout)
    tip.receive_sample_sheet(shell.supply_sample_sheet())
    print(f'{tip.module_name}: Removing samples lacking required files', file=sys.stdout)
    samples_cleared = tip.remove_invalid_samples(connect_from_module_name='core', taxonomy_only=True)
    tip.save_removed()
    if samples_cleared == 1:
        print(f'{tip.module_name}: No samples to process - proceeding to shape', file=sys.stdout)
        print(f'{shape.module_name}: Checking required input', file=sys.stdout)                                                     
        shape.receive_sample_sheet(shell.supply_sample_sheet())
        print(f'{shape.module_name}: Removing samples lacking required files', file=sys.stdout)
        samples_cleared = shape.remove_invalid_samples(connect_from_module_name='shell')

        # Running shape
        print(f'{shape.module_name}: Adding targets to config file', file=sys.stdout)
        shape.fill_input_dict(substring_list=None, mixed=True, empty=True)               #empty sample sheet due to filtering of invalid samples
        shape.fill_target_list(mixed=True, empty=True)
        if args.nanopore_only or args.fasta:
            print(f'{shape.module_name}: {"Nanopore" if args.nanopore_only else "Fasta"} mode requested - configuring', file=sys.stdout)
            shape.target_list = [t for t in shape.target_list if "fastp" not in t and "host_filtering" not in t]
            shape.snakefile_path = module_data['snakefiles']['shape_ont']
        shape.add_module_targets()
        shape.config_cluster()
        shape.write_module_config()
        print(f'{shape.module_name}: Launching snakemake', file=sys.stdout)
        shape.run_module(job_count=num_jobs)

        print(f'{shape.module_name}: Verifying output completion', file=sys.stdout)
        shape.check_module_output(mixed=True)
        shape.write_sample_sheet()
        if shape.pack_output:
            print('Packing_output')
            shape.output_path = tip.output_path
            shape.fold_output()
        shape.set_permissions()
        sys.exit(f'{shape.module_name}: Finished')
    else:
        # Running tip
        print(f'{tip.module_name}: Adding targets to config file', file=sys.stdout)
        tip.fill_input_dict(substring_list=None, pattern_path='ONT')
        tip.add_fasta_samples()
        tip.write_sample_sheet()
        tip.fill_target_list(taxonomy_based=True)
        if args.nanopore_only or args.fasta:
            print(f'{tip.module_name}: {"Nanopore" if args.nanopore_only else "Fasta"} mode requested - configuring', file=sys.stdout)
            tgt_exc_list = [
                "-predictResults.txt",
                "_SeqSero.tsv",
                "_stecfinder.tsv",
                "_seroba.tsv",
                "_lrefinder/lrefinder.tsv",
                "_shigatyper/shigatyper.tsv"
            ]
            tip.target_list = [t for t in tip.target_list if not any([tag in t for tag in tgt_exc_list])]
            tip.snakefile_path = module_data['snakefiles']['tip_ont']
        tip.add_module_targets()
        tip.config_cluster()
        tip.write_module_config()
        print(f'{tip.module_name}: Launching snakemake', file=sys.stdout)
        tip.run_module(job_count=num_jobs)

        print(f'{tip.module_name}: Verifying output completion', file=sys.stdout)
        tip.check_module_output()
        tip.write_sample_sheet()
        print(f'{tip.module_name}: Finished', file=sys.stdout)

    # Connecting tip & core to shape
    print(f'{shape.module_name}: Checking required input', file=sys.stdout)
    shape.receive_sample_sheet(tip.supply_sample_sheet())
    print(f'{shape.module_name}: Removing samples lacking required files', file=sys.stdout)
    samples_cleared = shape.remove_invalid_samples(connect_from_module_name='core')
    shape.removed_samples = tip.removed_samples
    if samples_cleared == 1: 
        if tip.pack_output: tip.fold_output()
        print(f'{shape.module_name}: Missing files requested by bact_shell', file=sys.stderr)
        sys.exit(1)

    # Running shape
    print(f'{shape.module_name}: Adding targets to config file', file=sys.stdout)
    shape.fill_input_dict(substring_list=None, mixed=True)
    shape.fill_target_list(mixed=True)
    if args.nanopore_only or args.fasta:
        print(f'{shape.module_name}: {"Nanopore" if args.nanopore_only else "Fasta"} mode requested - configuring', file=sys.stdout)
        tgt_exc_list = [
            "fastp",
            "host_filtering",
            "_SeqSero_std.csv",
            "_stecfinder_std.csv",
            "_seroba_std.csv",
            "-predictResults_std.csv"
        ]
        shape.target_list = [t for t in shape.target_list if not any([tag in t for tag in tgt_exc_list])]
        shape.snakefile_path = module_data['snakefiles']['shape_ont']
    shape.add_module_targets()
    shape.config_cluster()
    shape.write_module_config()
    print(f'{shape.module_name}: Launching snakemake', file=sys.stdout)
    shape.run_module(job_count=num_jobs)

    print(f'{shape.module_name}: Verifying output completion', file=sys.stdout)
    shape.check_module_output(mixed=True)
    shape.write_sample_sheet()
    if shape.pack_output:
        shape.output_path = tip.output_path
        tip.fold_output()
    shape.set_permissions()

    print(f'{shape.module_name}: Finished', file=sys.stdout)

    # Add sample_id and job name to log
    return shape.output_path


def run_merge(args, num_jobs):
    '''Wrapper function to combine outputs from multiple folders and run all modules on the result.'''
    if args.merge_from is None:
        print('Please specify 1 or more paths to --merge_from to run `merge` mode',file=sys.stderr)
        sys.exit(1)
    if args.output_dir is None:
        print('Please specify a path to --output_dir to run `merge` mode',file=sys.stderr)
        sys.exit(1)
    merge_inputs = args.merge_from
    merge_target = args.output_dir
    exclude_files = [
        'sample_sheet.csv', 
        'remove_samples_tip.csv', 
        'core_aggregated_taxonomy.json', 
        'cluster.yaml', 
        'config_cluste.yaml' , 
        'removed_samples_tip.csv'
        ]
    print(f'Merge: starting merging folder contents', file=sys.stdout)
    hk.merge_paths(src_list=merge_inputs, target_folder = merge_target, exclude_files = exclude_files)
    args.input = os.path.abspath(args.output_dir)
    if not args.input.endswith('/'):
        args.input = args.input+'/'
    args.skip_packing = False

    print(f'Merge: Launching the pipeline', file=sys.stdout)
    run_all(args, num_jobs)
    nl='\n'
    print(f'Succesfully merged\n\n{nl.join(merge_inputs)}\n\nin {merge_target}.', file=sys.stdout)
