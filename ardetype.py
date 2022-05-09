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
        print(f"\nStarting bact_core: read QC, host filtering, taxonomic classification and contig assembly.\n")
        run_core(args)
        assert os.path.isfile(f"{args.output_dir}/sample_sheet.csv"), f'Cannot find sample sheet in {args.output_dir} after bact_core stage is complete'
        print(f"\nStarting bact_shell:\n")
        run_shell(args)
        assert os.path.isfile(f"{args.output_dir}/sample_sheet.csv"), f'Cannot find sample sheet in {args.output_dir} after bact_shell stage is complete'
        # run_tip(args)
    elif args.mode == "core":
        run_core(args)
    elif args.mode == "shell":
        run_shell(args)
    elif args.mode == "tip":
        run_tip(args)
    if not args.skip_reporting:
        run_shape(args)