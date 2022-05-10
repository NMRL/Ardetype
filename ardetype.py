from ardetype_modules import run_core, run_shell, run_tip, run_shape
from ardetype_utilities import install_snakemake, parse_arguments
import os

"""
This is a wrapper script of ARDETYPE pipeline.
Date: 2022-05-09
Version: 0.1
"""

if __name__ == "__main__":
    args = parse_arguments()
    if args.install_snakemake:
        install_snakemake()
    if args.mode == "all":
        try:
            print(f"\nStarting bact_core: read QC, host filtering, taxonomic classification and contig assembly.\n")
            core_check_dict = run_core(args)
            print(f"\nStarting bact_shell:\n")
            shell_check_dict = run_shell(args)
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
    elif args.mode == "tip":
        run_tip(args)
    # if not args.skip_reporting:
    #     run_shape(args)