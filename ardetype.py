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
    if args.install_snakemake:
        hk.install_snakemake()

    if args.mode == "bucket":
        __config = hk.read_yaml("./config_files/yaml/config_modular.yaml")
        home_path = __config['home_dir']
        out_bucket_path = __config['output_bucket_path']
        del __config
        if args.input.endswith('/'):
            input_name = os.path.basename(os.path.dirname(args.input))
        else:
            input_name = os.path.basename(args.input)
        target_path = os.path.join(home_path, input_name)
        hk.move_folder(args.input, target_path)
        print(f'Moving from input bucket to {target_path}')
        args.input, args.output_dir = target_path, target_path

        output = run_all(args, num_jobs)

        if output.endswith('/'):
            target_name = os.path.basename(os.path.dirname(output))
        else:
            target_name = os.path.basename(output)

        target_path = os.path.join(home_path, target_name)
        target_bucket_path = os.path.join(out_bucket_path, target_name)
        print(f'Moving from {target_name} to {hk.move_folder(target_path, target_bucket_path)}')

    elif args.mode == "all":
        run_all(args, num_jobs)
    elif args.mode == 'log_analysis':
        hk.update_log_history(pipeline_name='ardetype')
        hk.update_log_summary(notebook_path='./subscripts/downstream/log_summary.ipynb', env_path='/mnt/home/groups/nmrl/cov_analysis/SARS-CoV2_assembly/tools/rbase_env/', output_dir='./')
    elif args.mode == 'merge':
        run_merge(args, num_jobs)
    else:
        sys.exit(f'Mode {args.mode} not supported, please use `--mode all`.')