#!/usr/bin/env python3

import subprocess as sp
import datetime
import os
import json
import argparse
from pathlib import Path
from subscripts.ardetype_utilities import Ardetype_housekeeper as hk

def generate_tool_versions(config_dict):
    """Replicates original parsing logic exactly."""

    def safe_run(func):
        try:
            return func()
        except Exception as e:
            return f"ERROR: {e}"

    tool_vers_map = {
            "plasmidfinder": 
                datetime.datetime.fromtimestamp(os.path.getmtime(config_dict["shell_tool_configs"]["plasmidfinder"]["plasmidfinder_sif"])).strftime('%Y-%m-%d'),
            "resfinder":
                sp.run(
                f'singularity --silent exec {config_dict["resfinder_sif"]} python -m resfinder --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip(),
            "virulencefinder":
                sp.run(f'singularity --silent exec {config_dict["shell_tool_configs"]["virulencefinder"]["virulencefinder_sif"]} python -m virulencefinder --version 2> /dev/null', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip(),
            "quast":
                sp.run(
                f'singularity --silent exec {config_dict["quast_sif"]} quast --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split('\n')[-1],
            "rgi":
                sp.run(
                f'singularity --silent exec {config_dict["rgi_sif"]} rgi main --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip(),
            "kraken2":
                sp.run(
                f'singularity --silent exec {config_dict["kraken2_sif"]} kraken2 --version',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split('\n')[0],
            "amrfinder+":
                sp.run(f'singularity --silent exec {config_dict["amrfinderplus_sif"]} amrfinder --version',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip(),
            "fastp": 
                sp.run(f'singularity --silent exec {config_dict["fastp_sif"]} fastp --version',
                stderr=sp.PIPE, shell=True).stderr.decode('utf-8').strip(),
            "mob-suite":
                sp.run(f'singularity --silent exec {config_dict["mob_suite_sif"]} mob_typer --version',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').split(" ")[-1].strip(),
            "mlst":
                sp.run(f'singularity --silent exec {config_dict["mlst_sif"]} mlst --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "shovill":
                sp.run(f'singularity --silent exec {config_dict["shovill_sif"]} shovill --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "meningotype":
                sp.run(f'singularity --silent exec {config_dict["meningotype_nmeningitidis_sif"]} meningotype --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "legsta":
                sp.run(f'singularity --silent exec {config_dict["legsta_lpneumophila_sif"]} legsta --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "hicap":
                sp.run(f'singularity --silent exec {config_dict["hicap_hinfluenzae_sif"]} hicap --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "kleborate":
                sp.run(f'singularity --silent exec {config_dict["kleborate_kpneumoniae_sif"]} kleborate --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "agrvate":
                sp.run(f'singularity --silent exec {config_dict["agrvate_saureus_sif"]} agrvate --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "spatyper":
                sp.run(f'singularity --silent exec {config_dict["spatyper_saureus_sif"]} spaTyper --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "staphopia-sccmec":
                sp.run(f'singularity --silent exec {config_dict["sccmec_saureus_sif"]} staphopia-sccmec --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "emmtyper":
                sp.run(f'singularity --silent exec {config_dict["emmtyper_spyogenes_sif"]} emmtyper --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "lissero":
                sp.run(f'singularity --silent exec {config_dict["lissero_lmonocytogenes_sif"]} lissero --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "sistr":
                sp.run(f'singularity --silent exec {config_dict["sistr_senterica_sif"]} sistr --version ',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "seqsero2":
                sp.run(f'singularity --silent exec {config_dict["seqsero2_senterica_sif"]} SeqSero2_package.py --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "ectyper":
                sp.run(f'singularity --silent exec {config_dict["ectyper_ecoli_sif"]} ectyper --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip(),
            "stecfinder":
                sp.run(f'singularity --silent exec {config_dict["stecfinder_ecoli_sif"]} stecfinder --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "seroba":
                sp.run(f'singularity --silent exec {config_dict["seroba_spneumoniae_sif"]} seroba version',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip(),
            "chewbbaca":
                sp.run(f'singularity --silent exec {config_dict["tip_tool_configs"]["chewbbaca"]["chewbbaca_sif"]} chewBBACA.py --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            'lrefinder':
                sp.run(f'singularity --silent exec  {config_dict["lrefinder_efaecium_efaecalis_sif"]} LRE-Finder.py -v', 
                stderr=sp.PIPE, shell=True).stderr.decode('utf-8').strip().split('\n')[-1],
            "shigatyper":
                sp.run(f'singularity --silent exec {config_dict["shigatyper_shigella_sif"]} shigatyper --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "flye":
                sp.run(f'singularity --silent exec {config_dict["flye_sif"]} flye --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "snikt":
                sp.run(f'singularity --silent exec {config_dict["snikt_sif"]} snikt.R --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "filtlong":
                sp.run(f'singularity --silent exec {config_dict["filtlong_sif"]} filtlong --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "circlator":
                sp.run(f'singularity --silent exec {config_dict["circlator_sif"]} circlator version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "polypolish":
                sp.run(f'singularity --silent exec {config_dict["polypolish_sif"]} polypolish --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1],
            "polca":"4.1.0",
            "medaka": sp.run(f'singularity --silent exec {config_dict["medaka_sif"]} medaka --version 2> /dev/null', stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split('\n')[-1].split()[1],
            "capybara":"1.0",
            "plasmidtype_abaumannii":"2022",
            "hybracter": "0.11.0"
        }

    return tool_vers_map


def main():
    parser = argparse.ArgumentParser(
        description="Generate static JSON with Singularity tool versions"
    )
    parser.add_argument(
        "-o", "--output",
        default="tool_versions.json",
        help="Output JSON file"
    )
    args = parser.parse_args()

    ardetype_path = str(Path(os.path.abspath('./')))
    config_dict = hk.read_yaml(os.path.join(ardetype_path,'config_files/yaml/config_modular_local.yaml'))

    versions = generate_tool_versions(config_dict)

    Path(args.output).write_text(json.dumps(versions, indent=2))
    print(f"Saved tool versions → {args.output}")


if __name__ == "__main__":
    main()
