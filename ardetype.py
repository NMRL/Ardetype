import sys, os
from pathlib import Path
sys.path.insert(0, f'{os.path.dirname(str(Path(__file__).absolute()))}/subscripts/')
from subscripts.ardetype_modules import run_all, run_merge
from subscripts.ardetype_utilities import Ardetype_housekeeper as hk

"""
This is a wrapper script of ARDETYPE pipeline.
Date: 2025-02-24
Version: 1.0.0
"""

if __name__ == "__main__":
    args = hk.parse_arguments(hk.read_json_dict('./config_files/json/argument_data.json'))
    num_jobs = args.num_jobs

    if args.mode == "all":
        run_all(args, num_jobs)
    elif args.mode == 'merge':
        run_merge(args, num_jobs)
    else:
        sys.exit(f'Mode {args.mode} not supported, please use `--mode all`.')