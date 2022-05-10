from os import remove
from ardetype_modules import run_core, run_shell, run_tip, run_shape
from ardetype_utilities import install_snakemake, parse_arguments, remove_invalid_samples

"""
This is a wrapper script of ARDETYPE pipeline.
Date: 2022-05-10
Version: 0.2
"""

if __name__ == "__main__":
    args = parse_arguments()
    if args.install_snakemake:
        install_snakemake()
    if args.mode == "all":
        try:
            print(f"\nStarting bact_core: read QC, host filtering, taxonomic classification and contig assembly.\n")
            core_sample_sheet = run_core(args)
            samples_cleared = remove_invalid_samples(core_sample_sheet,"shell",args.output_dir)
            print(f"\nStarting bact_shell:\n")
            assert samples_cleared is None, "\nNo samples can be processed by bact_shell module - required files are missing.\n"
            shell_sample_sheet = run_shell(args)
        except AssertionError as msg:
            print(f"{msg}")
        
        # run_tip(args)
    elif args.mode == "core":
        try:
            run_core(args)
        except AssertionError as msg:
            print(f"{msg}")
    elif args.mode == "shell":
        run_shell(args)
    # elif args.mode == "tip":
    #     run_tip(args)
    # if not args.skip_reporting:
    #     run_shape(args)