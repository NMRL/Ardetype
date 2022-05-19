from ardetype_modules import Module
from ardetype_utilities import install_snakemake, parse_arguments, read_json_dict
from pathlib import Path
import os

"""
This is a wrapper script of ARDETYPE pipeline.
Date: 2022-05-19
Version: 0.4
"""

if __name__ == "__main__":
    args = parse_arguments()
    args.output_dir = f"{os.path.abspath(args.output_dir)}/"
    num_jobs = args.num_jobs
    module_data = read_json_dict(f'{os.path.dirname(Path( __file__ ).absolute())}/module_data.json')

    if args.install_snakemake:
        install_snakemake()

    if args.mode == "all":
        #Running bact_core
        core = Module(
            module_name='core',
            input_path=args.input,
            module_config=args.config,
            output_path=args.output_dir,
            run_mode=args.module_jobs,
            job_name=module_data['core']['job_name'],
            patterns=module_data['core']['patterns'],
            targets=module_data['core']['targets'],
            requests=module_data['core']['requests'],
            snakefile_path=module_data['snakefiles']['core'],
            cluster_config_path=module_data['cluster_config']
            )
        core.fill_input_dict()
        core.fill_sample_sheet()
        core.make_output_dir()
        core.write_sample_sheet()
        core.fill_target_list()
        core.add_module_targets()
        core.add_output_dir()
        core.write_module_config()
        core.run_module(job_count=num_jobs)
        core.check_module_output()
        core.write_sample_sheet()

        #Running bact_shell
        shell = Module(
            module_name='shell', 
            input_path=core.output_path, 
            module_config=core.config_file, 
            output_path=args.output_dir, 
            run_mode=args.module_jobs,
            job_name=module_data['shell']['job_name'],
            patterns=module_data['shell']['patterns'],
            targets=module_data['shell']['targets'],
            requests=module_data['shell']['requests'],
            snakefile_path=module_data['snakefiles']['shell'],
            cluster_config_path=module_data['cluster_config']
            )
            
        #Connecting core to shell
        shell.receive_sample_sheet(core.supply_sample_sheet())
        samples_cleared = shell.remove_invalid_samples(connect_from_module_name='core')
        if samples_cleared == 1: raise Exception('Missing files requested by bact_shell.')

        shell.fill_input_dict()
        shell.add_fasta_samples()
        shell.write_sample_sheet()
        shell.fill_target_list()
        shell.add_module_targets()
        shell.write_module_config()
        shell.run_module(job_count=num_jobs)
        shell.check_module_output()
        shell.write_sample_sheet()


    elif args.mode == 'core':
        core = Module(
            module_name='core',
            input_path=args.input,
            module_config=args.config,
            output_path=args.output_dir,
            run_mode=args.module_jobs,
            job_name=module_data['core']['job_name'],
            patterns=module_data['core']['patterns'],
            targets=module_data['core']['targets'],
            requests=module_data['core']['requests'],
            snakefile_path=module_data['snakefiles']['core'],
            cluster_config_path=module_data['cluster_config']
        )
        core.fill_input_dict()
        core.fill_sample_sheet()
        core.make_output_dir()
        core.write_sample_sheet()
        core.fill_target_list()
        core.add_module_targets()
        core.add_output_dir()
        core.write_module_config()
        core.run_module(job_count=num_jobs)
        core.check_module_output()
        core.write_sample_sheet()


    elif args.mode == "shell":
        shell = Module(
            module_name='shell', 
            input_path=args.input,
            module_config=args.config, 
            output_path=args.output_dir, 
            run_mode=args.module_jobs,
            job_name=module_data['shell']['job_name'],
            patterns=module_data['shell']['patterns'],
            targets=module_data['shell']['targets'],
            requests=module_data['shell']['requests'],
            snakefile_path=module_data['snakefiles']['shell'],
            cluster_config_path=module_data['cluster_config']
        )
        shell.fill_input_dict()
        shell.fill_sample_sheet()
        shell.make_output_dir()
        shell.write_sample_sheet()
        shell.fill_target_list()
        shell.add_module_targets()
        shell.add_output_dir()
        shell.write_module_config()
        shell.run_module(job_count=num_jobs)
        shell.check_module_output()
        shell.write_sample_sheet()
