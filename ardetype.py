from modules import run_core
from ardetype_utilities import install_snakemake, parse_arguments

"""
This is a wrapper script of ARDETYPE pipeline.
Date: 2022-05-09
Version: 0.1
"""

if __name__ == "__main__":
    args = parse_arguments()
    if args.install_snakemake:
        install_snakemake()
    if args.mode == "core":
        run_core(args)
    else:
        print('Sorry, other options not supported yet.')
