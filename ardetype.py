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
    # Keep singularity image extraction and cache off the (small) root filesystem /tmp;
    # inherited by snakemake and the cluster jobs it submits. Pre-set values take precedence.
    pipeline_tmp = f'{os.path.dirname(str(Path(__file__).absolute()))}/tmp'
    for var, sub in [('SINGULARITY_TMPDIR', 'singularity_tmp'), ('APPTAINER_TMPDIR', 'singularity_tmp'),
                     ('SINGULARITY_CACHEDIR', 'singularity_cache'), ('APPTAINER_CACHEDIR', 'singularity_cache')]:
        os.environ.setdefault(var, f'{pipeline_tmp}/{sub}')
        os.makedirs(os.environ[var], exist_ok=True)

    # Put the singularity selector (subscripts/bin/singularity) first on PATH: it prefers a system-wide
    # singularity and falls back to the conda one; inherited by snakemake and the cluster jobs it submits.
    singularity_selector_dir = f'{os.path.dirname(str(Path(__file__).absolute()))}/subscripts/bin'
    if os.environ.get('PATH', '').split(':')[0] != singularity_selector_dir:
        os.environ['PATH'] = f"{singularity_selector_dir}:{os.environ.get('PATH', '')}"

    args = hk.parse_arguments(hk.read_json_dict('./config_files/json/argument_data.json'))
    num_jobs = args.num_jobs

    if args.run_mode == "all":
        run_all(args, num_jobs)
    elif args.run_mode == 'merge':
        run_merge(args, num_jobs)
    else:
        sys.exit(f'Mode {args.mode} not supported, please use `--mode all`.')