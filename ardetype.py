import sys, os
from pathlib import Path
sys.path.insert(0, f'{os.path.dirname(str(Path(__file__).absolute()))}/subscripts/')
from ardetype_modules import run_all, run_core, run_shell
from ardetype_utilities import Housekeeper as hk

"""
This is a wrapper script of ARDETYPE pipeline.
Date: 2022-05-30
Version: 1.0
"""

if __name__ == "__main__":
    args = hk.parse_arguments(hk.read_json_dict('./config_files/json/argument_data.json'))
    num_jobs = args.num_jobs
    if args.install_snakemake:
        hk.install_snakemake()
    if args.mode == "all":
        run_all(args, num_jobs)
    elif args.mode == 'core':
        run_core(args,num_jobs)
    elif args.mode == "shell":
        run_shell(args, num_jobs)