from subscripts.ardetype_utilities import Ardetype_housekeeper as hk
import os, warnings, sys
from pathlib import Path
sys.path.insert(0, os.path.dirname(Path(__file__).absolute()))

from src.modules import Module

#Suppressing pandas warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)


#Reading data used to build module objects
ardetype_path = os.path.dirname(Path(__file__).parents[0].absolute())
module_data   = hk.read_json_dict(f'{ardetype_path}/config_files/json/module_data.json')


###############################################################
# Extensing base Module class to fit the needs of the pipeline
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
            sys.exit(f'Input must be included to run the pipeline in any mode except log_analysis.')
        super(Ardetype_module, self).__init__(*args, **kwargs) #running method as it is defined in the base class


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
            dry_run             = args.dry_run,
            force_all           = args.force_all,
            rule_graph          = args.rule_graph,
            pack_output         = args.pack_output,
            unpack_output       = args.unpack_output,
            retry_times         = args.retry_times,
            job_name            = module_data['core']['job_name'],
            patterns            = module_data['core']['patterns'],
            targets             = module_data['core']['targets'],
            requests            = module_data['core']['requests'],
            snakefile_path      = module_data['snakefiles']['core'],
            cluster_config_path = module_data['cluster_config']
            )
    shell = Ardetype_module(
        module_name         = 'shell', 
        input_path          = core.output_path, 
        module_config       = core.config_file, 
        output_path         = args.output_dir, 
        run_mode            = args.submit_modules,
        dry_run             = args.dry_run,
        force_all           = args.force_all,
        rule_graph          = args.rule_graph,
        pack_output         = args.pack_output,
        unpack_output       = args.unpack_output,
        retry_times         = args.retry_times,
        job_name            = module_data['shell']['job_name'],
        patterns            = module_data['shell']['patterns'],
        targets             = module_data['shell']['targets'],
        requests            = module_data['shell']['requests'],
        snakefile_path      = module_data['snakefiles']['shell'],
        cluster_config_path = module_data['cluster_config']
        )
    tip = Ardetype_module(
        module_name         = 'tip', 
        input_path          = core.output_path,
        module_config       = shell.config_file, 
        output_path         = args.output_dir, 
        run_mode            = args.submit_modules,
        dry_run             = args.dry_run,
        force_all           = args.force_all,
        rule_graph          = args.rule_graph,
        pack_output         = args.pack_output,
        unpack_output       = args.unpack_output,
        retry_times         = args.retry_times,
        job_name            = module_data['tip']['job_name'],
        patterns            = module_data['tip']['patterns'],
        targets             = module_data['tip']['targets'],
        requests            = module_data['tip']['requests'],
        snakefile_path      = module_data['snakefiles']['tip'],
        cluster_config_path = module_data['cluster_config']
    )
    shape = Ardetype_module(
        module_name         = 'shape', 
        input_path          = core.output_path,
        module_config       = tip.config_file, 
        output_path         = args.output_dir, 
        run_mode            = args.submit_modules,
        dry_run             = args.dry_run,
        force_all           = args.force_all,
        rule_graph          = args.rule_graph,
        pack_output         = args.pack_output,
        unpack_output       = args.unpack_output,
        retry_times         = args.retry_times,
        job_name            = module_data['shape']['job_name'],
        patterns            = module_data['shape']['patterns'],
        targets             = module_data['shape']['targets'],
        requests            = module_data['shape']['requests'],
        snakefile_path      = module_data['snakefiles']['shape'],
        cluster_config_path = module_data['cluster_config']
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
    core.files_to_wd(redirect_filter={"001.fastq.gz":core.output_path})
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

    # Connecting shell & core to tip/shape
    tip.receive_sample_sheet(shell.supply_sample_sheet())
    samples_cleared = tip.remove_invalid_samples(connect_from_module_name='core')
    tip.save_removed()
    if samples_cleared == 1:                                                            
        shape.receive_sample_sheet(shell.supply_sample_sheet())
        samples_cleared = shape.remove_invalid_samples(connect_from_module_name='shell')

        # Running shape
        shape.fill_input_dict(substring_list=None, mixed=True, empty=True)               #empty sample sheet due to filtering of invalid samples
        shape.fill_target_list(mixed=True, empty=True)
        shape.add_module_targets()
        shape.write_module_config()
        try:
            shape.run_module(job_count=num_jobs)
        except Exception as e:
            raise e
        shape.check_module_output(mixed=True)
        shape.write_sample_sheet()
        if shape.pack_output: tip.fold_output()
        shape.set_permissions()
        sys.exit("bact_shape finished")
    else:
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

    # Connecting tip & core to shape
    shape.receive_sample_sheet(tip.supply_sample_sheet())
    samples_cleared = shape.remove_invalid_samples(connect_from_module_name='core')
    shape.removed_samples = tip.removed_samples
    if samples_cleared == 1: 
        if tip.pack_output: tip.fold_output()
        raise Exception('Missing files requested by bact_shape.')

    # Running shape
    shape.fill_input_dict(substring_list=None, mixed=True)
    shape.fill_target_list(mixed=True)
    shape.add_module_targets()
    shape.write_module_config()
    try:
        shape.run_module(job_count=num_jobs)
    except Exception as e:
        raise e
    shape.check_module_output(mixed=True)
    shape.write_sample_sheet()
    if shape.pack_output: tip.fold_output()
    shape.set_permissions()

    #Housekeeping the log files
    hk.asign_perm_rec(f"{ardetype_path}/ardetype_job_logs/")
    hk.name_job_logs('ardetype')
    if args.clean_job_logs:
        hk.remove_old_files(f"{ardetype_path}/ardetype_job_logs/")



def run_core(args, num_jobs):
    '''Wrapper function to run only core module.'''
    core = Ardetype_module(
        module_name         = 'core',
        input_path          = args.input,
        module_config       = args.config,
        output_path         = args.output_dir,
        run_mode            = args.submit_modules,
        dry_run             = args.dry_run,
        force_all           = args.force_all,
        rule_graph          = args.rule_graph,
        pack_output         = args.pack_output,
        unpack_output       = args.unpack_output,
        retry_times         = args.retry_times,
        job_name            = module_data['core']['job_name'],
        patterns            = module_data['core']['patterns'],
        targets             = module_data['core']['targets'],
        requests            = module_data['core']['requests'],
        snakefile_path      = module_data['snakefiles']['core'],
        cluster_config_path = module_data['cluster_config']
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
    core.set_permissions()

    #Housekeeping the log files
    hk.asign_perm_rec(f"{ardetype_path}/ardetype_job_logs/")
    hk.name_job_logs('ardetype')
    if args.clean_job_logs:
        hk.remove_old_files(f"{ardetype_path}/ardetype_job_logs/")


def run_shell(args, num_jobs):
    '''Wrapper function to run only shell module.'''
    shell = Ardetype_module(
        module_name         = 'shell', 
        input_path          = args.input,
        module_config       = args.config, 
        output_path         = args.output_dir, 
        run_mode            = args.submit_modules,
        dry_run             = args.dry_run,
        force_all           = args.force_all,
        rule_graph          = args.rule_graph,
        pack_output         = args.pack_output,
        unpack_output       = args.unpack_output,
        retry_times         = args.retry_times,
        job_name            = module_data['shell']['job_name'],
        patterns            = module_data['shell']['patterns'],
        targets             = module_data['shell']['targets'],
        requests            = module_data['shell']['requests'],
        snakefile_path      = module_data['snakefiles']['shell'],
        cluster_config_path = module_data['cluster_config']
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
    shell.set_permissions()

    #Housekeeping the log files
    hk.asign_perm_rec(f"{ardetype_path}/ardetype_job_logs/")
    hk.name_job_logs('ardetype')
    if args.clean_job_logs:
        hk.remove_old_files(f"{ardetype_path}/ardetype_job_logs/")