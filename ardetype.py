from ardetype_modules import run_all, run_core, run_shell
from ardetype_utilities import install_snakemake, parse_arguments

"""
This is a wrapper script of ARDETYPE pipeline.
Date: 2022-05-27
Version: 0.8
"""

if __name__ == "__main__":
    args = parse_arguments()
    num_jobs = args.num_jobs
    if args.install_snakemake:
        install_snakemake()
    if args.mode == "all":
        run_all(args, num_jobs)
    elif args.mode == 'core':
        run_core(args,num_jobs)
    elif args.mode == "shell":
        run_shell(args, num_jobs)