import sys, os, shutil
import yaml
from pathlib import Path
sys.path.insert(0, f'{os.path.dirname(str(Path(__file__).absolute()))}/subscripts/')
from subscripts.ardetype_modules import run_all, run_merge
from subscripts.ardetype_utilities import Ardetype_housekeeper as hk

"""
This is a wrapper script of ARDETYPE pipeline.
Date: 2026-11-16
Version: 1.0.1
"""

if __name__ == "__main__":
    # Keep singularity image extraction and cache off the (small) root filesystem /tmp;
    # inherited by snakemake and the cluster jobs it submits. Pre-set values take precedence.
    pipeline_tmp = f'{os.path.dirname(str(Path(__file__).absolute()))}/tmp'
    for var, sub in [('SINGULARITY_TMPDIR', 'singularity_tmp'), ('APPTAINER_TMPDIR', 'singularity_tmp'),
                     ('SINGULARITY_CACHEDIR', 'singularity_cache'), ('APPTAINER_CACHEDIR', 'singularity_cache')]:
        os.environ.setdefault(var, f'{pipeline_tmp}/{sub}')
        os.makedirs(os.environ[var], exist_ok=True)

    args = hk.parse_arguments(hk.read_json_dict('./config_files/json/argument_data.json'))

    # Singularity executable: taken from `singularity_bin` in the config file (see README, Configuration).
    # A link named `singularity` to it is put first on PATH, so snakemake and the cluster jobs it submits
    # (which inherit PATH) all use it. If `singularity_bin` is empty, the first `singularity` on PATH is used.
    with open(args.config) as fh:
        singularity_bin = (yaml.safe_load(fh) or {}).get('singularity_bin') or ''
    if singularity_bin:
        singularity_bin = os.path.abspath(os.path.expanduser(str(singularity_bin)))
        if not (os.path.isfile(singularity_bin) and os.access(singularity_bin, os.X_OK)):
            sys.exit(f"singularity_bin in {args.config} is not an executable file: {singularity_bin}")
        link_dir = f'{pipeline_tmp}/singularity_bin'
        os.makedirs(link_dir, exist_ok=True)
        tmp_link = f'{link_dir}/singularity.{os.getpid()}'
        os.symlink(singularity_bin, tmp_link)
        os.replace(tmp_link, f'{link_dir}/singularity')
        os.environ['PATH'] = f"{link_dir}:{os.environ.get('PATH', '')}"
    resolved = shutil.which('singularity')
    if resolved is None:
        sys.exit("No singularity executable found: set singularity_bin in the config file (see README, Configuration).")
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if not singularity_bin and conda_prefix and os.path.realpath(resolved).startswith(os.path.realpath(conda_prefix) + os.sep):
        print(f"WARNING: using the conda-provided singularity ({resolved}). It unpacks every image into "
              f"SINGULARITY_TMPDIR and is much slower than a system installation; set singularity_bin "
              f"in {args.config} if one is available (see README, Configuration).", file=sys.stderr)
    print(f"Using singularity: {os.path.realpath(resolved)}", file=sys.stderr)
    num_jobs = args.num_jobs

    if args.run_mode == "all":
        run_all(args, num_jobs)
    elif args.run_mode == 'merge':
        run_merge(args, num_jobs)
    else:
        sys.exit(f'Mode {args.mode} not supported, please use `--mode all`.')