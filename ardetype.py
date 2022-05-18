from ardetype_modules import Module
from ardetype_utilities import install_snakemake, parse_arguments, remove_invalid_samples
import os, pandas as pd

"""
This is a wrapper script of ARDETYPE pipeline.
Date: 2022-05-18
Version: 0.3
"""

if __name__ == "__main__":
    args = parse_arguments()
    args.output_dir = f"{os.path.abspath(args.output_dir)}/"

    if args.install_snakemake:
        install_snakemake()

    if args.mode == "all":
        #Running bact_core
        core = Module(module_name='core', input_path=args.input, module_config=args.config, output_path=args.output_dir, run_mode=args.module_jobs)
        core.fill_input_dict()
        core.fill_sample_sheet()
        core.make_output_dir()
        core.write_sample_sheet()
        core.fill_target_list()
        core.add_module_targets()
        core.add_output_dir()
        core.write_module_config()
        core.run_module()
        core.check_module_output()
        core.write_sample_sheet()

        #Running bact_shell
        shell = Module(module_name='shell', input_path=core.output_path, module_config=core.config_file, output_path=args.output_dir, run_mode=args.module_jobs)
        shell.fill_input_dict()
        samples_cleared = remove_invalid_samples(core.supply_sample_sheet(),"shell",args.output_dir)
        if not isinstance(samples_cleared, pd.DataFrame): raise Exception('Missing files, requested by bact_shell.')
        shell.receive_sample_sheet(samples_cleared)
        shell.add_fasta_samples()
        shell.write_sample_sheet()
        shell.fill_target_list()
        shell.add_module_targets()
        shell.write_module_config()
        shell.run_module()
        shell.check_module_output()
        shell.write_sample_sheet()


    elif args.mode == 'core':
        core = Module(module_name='core', input_path=args.input, module_config=args.config, output_path=args.output_dir, run_mode=args.module_jobs)
        core.fill_input_dict()
        core.fill_sample_sheet()
        core.make_output_dir()
        core.write_sample_sheet()
        core.fill_target_list()
        core.add_module_targets()
        core.add_output_dir()
        core.write_module_config()
        core.run_module()
        core.check_module_output()
        core.write_sample_sheet()


    elif args.mode == "shell":
        shell = Module(module_name='shell', input_path=args.input, module_config=args.config, output_path=args.output_dir, run_mode=args.module_jobs)
        shell.fill_input_dict()
        shell.fill_sample_sheet()
        shell.make_output_dir()
        shell.write_sample_sheet()
        shell.fill_target_list()
        shell.add_module_targets()
        shell.add_output_dir()
        shell.write_module_config()
        shell.run_module()
        shell.check_module_output()
        shell.write_sample_sheet()

    # elif args.mode == "tip":
    #     tip = Module(args, 'report')

    # if not args.skip_reporting:
    #   report = Module(args, 'report')