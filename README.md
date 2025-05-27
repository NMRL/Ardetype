# ARDETYPE

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/release/NMRL/Ardetype?color=Green)

> **ARDETYPE** is a modular, scalable pipeline for bacterial whole genome sequencing (WGS) analysis, supporting both **Illumina (paired-end)** and **Oxford Nanopore** reads. It is designed for clinical microbiology, public health, and research settings.

---

## Table of Contents

- [Features](#features)
- [Pipeline Modules](#pipeline-modules)
- [Supported Tools](#supported-tools)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)

---

## Features

- Supports both Illumina and Oxford Nanopore sequencing data
- Species-agnostic and species-specific analysis
- Modular, Snakemake-based workflow
- Automated quality control, assembly, classification, AMR profiling, plasmid typing, and reporting

---

## Pipeline Modules

ARDETYPE consists of four main modules:

| Module      | Description                                                                      |
|-------------|----------------------------------------------------------------------------------|
| `bact_core` | Core tasks: QC, contamination filtering, species ID, assembly                    |
| `bact_shell`| General tools: resistance gene & plasmid detection, assembly QC                  |
| `bact_tip`  | Species-specific typing and subtyping                                            |
| `shape`     | Aggregation, summarization, and report generation                                |

---

## Supported Tools

| Tool                                                                 | Module     |
|----------------------------------------------------------------------|------------|
| [agrvate](https://github.com/VishnuRaghuram94/AgrVATE)               | bact_tip   |
| [amrfinder+](https://github.com/ncbi/amr)                            | bact_shell |
| [capybara](https://github.com/Zhou-lab-SUDA/CAPYBARA)                | bact_tip   |
| [chewbbaca](https://github.com/B-UMMI/chewBBACA)                     | bact_tip   |
| [circlator](https://github.com/sanger-pathogens/circlator)           | bact_core  |
| [ectyper](https://github.com/phac-nml/ecoli_serotyping)              | bact_tip   |
| [emmtyper](https://github.com/MDU-PHL/emmtyper)                      | bact_tip   |
| [fastp](https://github.com/OpenGene/fastp)                           | bact_core  |
| [filtlong](https://github.com/rrwick/Filtlong)                       | bact_core  |
| [flye](https://github.com/fenderglass/Flye)                          | bact_core  |
| [hicap](https://github.com/scwatts/hicap)                            | bact_tip   |
| [kleborate](https://github.com/katholt/Kleborate)                    | bact_tip   |
| [kraken2](https://github.com/DerrickWood/kraken2)                    | bact_core  |
| [legsta](https://github.com/tseemann/legsta)                         | bact_tip   |
| [lissero](https://github.com/MDU-PHL/LisSero)                        | bact_tip   |
| [lrefinder](https://bitbucket.org/genomicepidemiology/lre-finder)    | bact_tip   |
| [medaka](https://github.com/nanoporetech/medaka)                     | bact_core  |
| [meningotype](https://github.com/MDU-PHL/meningotype)                | bact_tip   |
| [mlst](https://github.com/tseemann/mlst)                             | bact_shell |
| [mob-suite](https://github.com/phac-nml/mob-suite)                   | bact_shell |
| [plasmidfinder](https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/) | bact_shell |
| [polca](https://github.com/alekseyzimin/masurca)                     | bact_core  |
| [polypolish](https://github.com/rrwick/Polypolish)                   | bact_core  |
| [quast](https://github.com/ablab/quast)                              | bact_shell |
| [resfinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) | bact_shell |
| [rgi](https://github.com/arpcard/rgi)                                | bact_shell |
| [seqsero2](https://github.com/denglab/SeqSero2)                      | bact_tip   |
| [seroba](https://github.com/sanger-pathogens/seroba)                 | bact_tip   |
| [shigatyper](https://github.com/CFSAN-Biostatistics/shigatyper)      | bact_tip   |
| [shovill](https://github.com/tseemann/shovill)                       | bact_core  |
| [sistr](https://github.com/phac-nml/sistr_cmd)                       | bact_tip   |
| [snikt](https://github.com/piyuranjan/SNIKT)                         | bact_core  |
| [spatyper](https://github.com/HCGB-IGTP/spaTyper)                    | bact_tip   |
| [staphopia-sccmec](https://github.com/staphopia/staphopia-sccmec)    | bact_tip   |
| [stecfinder](https://github.com/LanLab/STECFinder)                   | bact_tip   |
| [virulencefinder](https://bitbucket.org/genomicepidemiology/virulencefinder/src/master/) | bact_shell |

---

## Installation

```bash
git clone -b standalone_rework https://github.com/NMRL/Ardetype
cd Ardetype
conda env create -f ardetype.yaml
```

**Download pre-packaged resources (databases & singularity images):**
```bash
rm -rf resources
wget <archive_url>
tar -xvzf resources.tar.gz
```

---

## Configuration

To enable **[rMLST](https://pubmlst.org/species-id)** support, configure API access to the PubMLST REST service by creating a file at `config_files/.env` with your credentials:

```env
CONSUMER_KEY=""
CONSUMER_SECRET=""
TEST_REST_URL="https://rest.pubmlst.org/db/pubmlst_rmlst_seqdef_kiosk/schemes/1"
SESSION_REST_URL="https://rest.pubmlst.org/db/pubmlst_test_seqdef"
TEST_WEB_URL="https://pubmlst.org/bigsdb?db=pubmlst_test_seqdef"
ACCESS_TOKEN_PATH="<path_for_access_token>"
```
- Request an API key pair from PubMLST and fill in `CONSUMER_KEY` and `CONSUMER_SECRET`.
- See [official documentation](https://bigsdb.readthedocs.io/en/latest/rest.html#authentication) for details.
- If rMLST is not configured, species will be inferred using Kraken2 with the MiniKraken database.

---

## Usage

**Typical Examples:**
```bash
# Illumina data
python ardetype.py -i <batch_path> -m all

# Nanopore data
python ardetype.py -i <batch_path> -m all -ont
```
- It is recommended to use the `screen` utility to run the pipeline in the background, especially for long-running jobs or remote sessions.

ARDETYPE integrates open-source tools for:

- Quality control, assembly, and taxonomic classification
- Antimicrobial resistance (AMR) prediction
- Plasmid reconstruction and typing
- Structured, integrated reporting

**Pipeline modules:**
- `bact_core`: QC, host read removal, assembly, classification
- `bact_shell`: assembly QC, AMR profiling, plasmid typing
- `bact_tip`: species-specific subtyping
- `bact_shape`: results summarization

**Command-line options:**

<details>
<summary>Click to expand</summary>

```text
optional arguments:
  -h, --help            Show this help message and exit
  -i, --input           Path to directory containing files to be analysed.
  -c, --config          Path to custom pipeline config file (YAML).
  -o, --output_dir      Directory for pipeline outputs.
  -r, --retry_times     Number of times to rerun failed jobs (default: 3).
  -n, --num_jobs        Max parallel Snakemake jobs (default: 12).
  -hpc, --hpc           Enable HPC execution mode.
  -sm, --submit_modules Run each module as separate HPC job (with --hpc).
  -f, --force_all       Force rerun of all rules.
  -sp, --skip_packing   Disable output packing.
  -ont, --nanopore_only Run exclusively on Nanopore reads.
  -mf, --merge_from     Merge outputs from multiple runs.
  -fr, --force_rules    Specify rules to force rerun.

required arguments:
  -m, --run_mode
    all   - Full pipeline (all modules)
    merge - Combine outputs from previous ARDETYPE runs
```
</details>