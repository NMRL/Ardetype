#!/usr/bin/env python3

import subprocess as sp
import datetime
import os
import re
import json
import argparse
import sys
from pathlib import Path

# Script lives in subscripts/; repo root is one level up
ardetype_path = str(Path(__file__).resolve().parents[1])
sys.path.insert(0, ardetype_path)
from subscripts.ardetype_utilities import Ardetype_housekeeper as hk

# Use the singularity selector (prefers system singularity over the slow conda build), also when run standalone
os.environ['PATH'] = f"{os.path.join(ardetype_path, 'subscripts', 'bin')}:{os.environ.get('PATH', '')}"

def safe_run(func):
    '''Runs a version query; a failure or empty output yields an "ERROR: ..." value instead of aborting the script.'''
    try:
        value = func()
    except Exception as e:
        value = f"ERROR: {type(e).__name__}: {e}"
    else:
        if not str(value).strip():
            value = "ERROR: empty version output"
    if str(value).startswith("ERROR"):
        print(f"WARNING: version query failed -> {value}", file=sys.stderr)
    return value


def query_version(sif:str, cmd:str, pattern:str) -> str:
    '''Runs cmd inside sif and returns the first regex group found in combined stdout+stderr.'''
    out = sp.run(f'singularity --silent exec {sif} {cmd}',
                 stdout=sp.PIPE, stderr=sp.STDOUT, shell=True).stdout.decode('utf-8')
    match = re.search(pattern, out, re.MULTILINE)
    if match is None:
        raise ValueError(f"version pattern not found in output of '{cmd}'")
    return match.group(1)


def image_build_date(sif:str) -> str:
    '''Returns the build date (YYYY-MM-DD) from the image label org.label-schema.build-date,
    e.g. "Tuesday_20_May_2025_13:13:28_EEST" -> "2025-05-20".'''
    out = sp.run(f'singularity inspect --labels {sif}', stdout=sp.PIPE, stderr=sp.STDOUT, shell=True).stdout.decode('utf-8')
    match = re.search(r'org\.label-schema\.build-date:\s*\w+_(\d{1,2}_\w+_\d{4})', out)
    if match is None:
        raise ValueError(f"build-date label not found for {sif}")
    return datetime.datetime.strptime(match.group(1), '%d_%B_%Y').strftime('%Y-%m-%d')


def generate_image_db_versions(config_dict):
    """Versions of databases bundled inside images (no external DB path is passed to these tools)."""
    return {
        # AMRFinderPlus runs without -d -> uses the database shipped in the image
        "amrfinder+": safe_run(lambda: query_version(config_dict["amrfinderplus_sif"], 'amrfinder -V', r'Database version:\s*(\S+)')),
        # RGI runs without --local -> uses the CARD data bundled with the rgi package (card.json "_version")
        "rgi": safe_run(lambda: query_version(config_dict["rgi_sif"],
            '''python -c "import app, os, json; print('CARD', json.load(open(os.path.join(os.path.dirname(app.__file__), '_data', 'card.json')))['_version'])"''',
            r'^CARD\s+(\S+)')),
    }


def generate_tool_versions(config_dict):
    """Replicates original parsing logic exactly."""

    # CAPYBARA has a single release tag (Capybara.1.0, 2024-04-12, commit 6f2146e); images are built from main afterwards.
    # The image has no git metadata, so the upstream commit is identified by the md5 of capy.py
    # (https://github.com/Zhou-lab-SUDA/CAPYBARA/commits/main). Each md5 maps to the commit that introduced that capy.py;
    # extend this table when a newer image is built.
    CAPYBARA_RELEASE = "1.0"
    CAPYBARA_CAPY_MD5 = {
        "a86ca74a3b0243f501312a0a5416b510": "119c20b",  # 2024-07-17
        "22ea13ab34427e6a6d23e6c752788783": "c9ba7b1",  # 2025-09-28
        "f286c4346ce2c0304c32830b118a4d7d": "3daa63c",  # 2026-03-25
        "4adb10104c4a96e615bdb15157f32e72": "bf4af18",  # 2026-08-27
    }

    def capybara_version() -> str:
        md5 = query_version(config_dict["tip_tool_configs"]["capybara_abaumannii"]["capybara_sif"],
                            'md5sum /capybara/capy.py', r'^([0-9a-f]{32})')
        commit = CAPYBARA_CAPY_MD5.get(md5, f"unknown, capy.py md5 {md5}")
        return f"{CAPYBARA_RELEASE} (commit {commit})"

    def download_date(path:str) -> str:
        '''Returns the modification date (= download date) of a database file as YYYY-MM-DD.'''
        return datetime.datetime.fromtimestamp(os.path.getmtime(path)).strftime('%Y-%m-%d')

    def abp_version() -> str:
        # AcinetobacterPlasmidTyping is a sequence database without tags/releases -> report download date
        # key is 'database' (not '*_database'), so join_sif_paths leaves it relative -> resolve against repo root
        return download_date(os.path.join(ardetype_path, config_dict["tip_tool_configs"]["abaumanniiplasmidtype"]["database"]))

    tool_vers_map = {
            # staphb/plasmidfinder:latest predates upstream versioning -> report the image build date
            "plasmidfinder":
                safe_run(lambda: image_build_date(config_dict["shell_tool_configs"]["plasmidfinder"]["plasmidfinder_sif"])),
            "resfinder":
                safe_run(lambda: sp.run(
                f'singularity --silent exec {config_dict["resfinder_sif"]} python -m resfinder --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip()),
            "virulencefinder":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["shell_tool_configs"]["virulencefinder"]["virulencefinder_sif"]} python -m virulencefinder --version 2> /dev/null', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip()),
            "quast":
                safe_run(lambda: sp.run(
                f'singularity --silent exec {config_dict["quast_sif"]} quast --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split('\n')[-1]),
            "rgi":
                safe_run(lambda: sp.run(
                f'singularity --silent exec {config_dict["rgi_sif"]} rgi main --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip()),
            "kraken2":
                safe_run(lambda: sp.run(
                f'singularity --silent exec {config_dict["kraken2_sif"]} kraken2 --version',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split('\n')[0]),
            "amrfinder+":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["amrfinderplus_sif"]} amrfinder --version',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip()),
            "fastp": 
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["fastp_sif"]} fastp --version',
                stderr=sp.PIPE, shell=True).stderr.decode('utf-8').strip()),
            "mob-suite":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["mob_suite_sif"]} mob_typer --version',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').split(" ")[-1].strip()),
            "mlst":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["mlst_sif"]} mlst --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            "shovill":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["shovill_sif"]} shovill --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            "meningotype":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["meningotype_nmeningitidis_sif"]} meningotype --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            # L. pneumophila SBT
            "el_gato": safe_run(lambda: query_version(config_dict["elgato_legionella_sif"], 'el_gato.py --version', r'el_gato version:\s*(\S+)')),
            "hicap":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["hicap_hinfluenzae_sif"]} hicap --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            "kleborate":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["kleborate_kpneumoniae_sif"]} kleborate --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            "agrvate":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["agrvate_saureus_sif"]} agrvate --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            # spaTyper --version reports a stale hard-coded 0.2.1; the installed package metadata is correct
            "spatyper": safe_run(lambda: query_version(config_dict["spatyper_saureus_sif"], 'pip show spatyper', r'^Version:\s*(\S+)')),
            "staphopia-sccmec":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["sccmec_saureus_sif"]} staphopia-sccmec --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            "emmtyper":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["emmtyper_spyogenes_sif"]} emmtyper --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            "lissero":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["lissero_lmonocytogenes_sif"]} lissero --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            "sistr":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["sistr_senterica_sif"]} sistr --version ',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            "seqsero2":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["seqsero2_senterica_sif"]} SeqSero2_package.py --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            "ectyper":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["ectyper_ecoli_sif"]} ectyper --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip()),
            "stecfinder":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["stecfinder_ecoli_sif"]} stecfinder --version 2> /dev/null',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            "seroba":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["seroba_spneumoniae_sif"]} seroba version',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip()),
            "chewbbaca":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["tip_tool_configs"]["chewbbaca"]["chewbbaca_sif"]} chewBBACA.py --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            'lrefinder':
                safe_run(lambda: sp.run(f'singularity --silent exec  {config_dict["lrefinder_efaecium_efaecalis_sif"]} LRE-Finder.py -v', 
                stderr=sp.PIPE, shell=True).stderr.decode('utf-8').strip().split('\n')[-1]),
            "shigatyper":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["shigatyper_shigella_sif"]} shigatyper --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            "flye":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["flye_sif"]} flye --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            "snikt":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["snikt_sif"]} snikt.R --version',
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            "fastplong": safe_run(lambda: query_version(config_dict["fastplong_sif"], 'fastplong --version', r'fastplong\s+(\S+)')),
            "filtlong":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["filtlong_sif"]} filtlong --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            "circlator":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["circlator_sif"]} circlator version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            "polypolish":
                safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["polypolish_sif"]} polypolish --version', 
                stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split(' ')[-1]),
            # POLCA is currently not used -> kept static
            "polca":"4.1.0",
            "medaka": safe_run(lambda: sp.run(f'singularity --silent exec {config_dict["medaka_sif"]} medaka --version 2> /dev/null', stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip().split('\n')[-1].split()[1]),
            # capy.py has no version option -> identify upstream commit by checksum
            "capybara": safe_run(capybara_version),
            # Acinetobacter rep typing database -> download date of the database file
            "plasmidtype_abaumannii": safe_run(abp_version),
            # hybracter has no working --version flag; read the installed package metadata
            "hybracter": safe_run(lambda: query_version(config_dict["hybracter_sif"], 'pip show hybracter', r'^Version:\s*(\S+)')),
            # KronaTools has no --version flag; the banner printed without arguments carries the version
            "krona": safe_run(lambda: query_version(config_dict["krona_sif"], 'ktImportText', r'KronaTools\s+(\S+)')),
            "ngmaster": safe_run(lambda: query_version(config_dict["tip_tool_configs"]["ngmaster_neisseria"]["ngmaster_sif"], 'ngmaster --version', r'ngmaster\s+(\S+)')),
            "hamronization": safe_run(lambda: query_version(config_dict["hamronization_sif"], 'hamronize --version', r'hamronize\s+(\S+)')),
            "cgmlst-dists": safe_run(lambda: query_version(config_dict["cgmlst_dists_sif"], 'cgmlst-dists -v', r'cgmlst-dists\s+(\S+)')),
            # BLAST+ used for A. baumannii plasmid rep typing
            "blast": safe_run(lambda: query_version(config_dict["tip_tool_configs"]["abaumanniiplasmidtype"]["abp_sif"], 'blastn -version', r'blastn:\s*(\S+)'))
        }

    return tool_vers_map


def main():
    parser = argparse.ArgumentParser(
        description="Generate static JSON with Singularity tool versions"
    )
    parser.add_argument(
        "-o", "--output",
        default=os.path.join(ardetype_path, "config_files/json/tool_versions.json"),
        help="Output JSON file (default: config_files/json/tool_versions.json)"
    )
    parser.add_argument(
        "-d", "--db-output",
        default=os.path.join(ardetype_path, "config_files/json/image_db_versions.json"),
        help="Output JSON file for databases bundled in images (default: config_files/json/image_db_versions.json)"
    )
    args = parser.parse_args()

    config_dict = hk.read_yaml(os.path.join(ardetype_path,'config_files/yaml/config_modular_local.yaml'))
    # Resolve relative image/db paths against the repo root, so the script works from any cwd
    for substr in ["_database", "_db", "_sif"]:
        hk.join_sif_paths(config_dict, substr, ardetype_path)

    versions = generate_tool_versions(config_dict)

    Path(args.output).write_text(json.dumps(versions, indent=2))
    print(f"Saved tool versions → {args.output}")

    db_versions = generate_image_db_versions(config_dict)
    Path(args.db_output).write_text(json.dumps(db_versions, indent=2))
    print(f"Saved image database versions → {args.db_output}")


if __name__ == "__main__":
    main()
