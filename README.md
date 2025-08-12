# ARDETYPE

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/release/NMRL/Ardetype?color=Green)

> **ARDETYPE** is a modular, scalable pipeline for bacterial whole-genome sequencing (WGS) analysis. It supports **Illumina (paired-end)** and **Oxford Nanopore** data in `fastq.gz` format, as well as **PacBio** reads and pre-assembled genomes in FASTA format.<br>Designed for clinical microbiology and public health applications, it enables streamlined and reproducible analysis of bacterial isolates through both species-agnostic and species-specific workflows.


---

## Features

- Species-agnostic and species-specific analysis of bacterial WGS data
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
| `bact_shape`| Aggregation, summarization, and report generation                                |

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
| [hybracter](https://github.com/gbouras13/hybracter) | bact_core |
| [cgmlst-dists](https://github.com/tseemann/cgmlst-dists) | bact_shape |
| [AcinetobacterPlasmidTyping](https://github.com/MehradHamidian/AcinetobacterPlasmidTyping) | bact_tip |
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
- If rMLST `is not configured`, species will be inferred using Kraken2 with the MiniKraken database.
- If rMLST `is configured`, species will be determined by a voting system that gives precedence to rMLST over Kraken2.

The `config_files/yaml/config_modular_local.yaml` file can be customized—for example, to define a custom location for the scratch directory used for intensive I/O operations.

- Relative paths should be specified with respect to the root of the `Ardetype` repository folder.
- Absolute paths are also supported, as long as they are accessible from the machine running the pipeline.

---

## Usage

**Input Structure Examples:**

- `isolate` stands for sample ID.
- Having `_` as part of sample ID may cause unexpected behavior.
- `isolate` in total `must shorter than 27 characters`, otherwise chewbbaca (cgMLST) will fail.
```bash
illumina/
├── isolate_R1_001.fastq.gz
└── isolate_R2_001.fastq.gz

ont/
└── isolate_ONT.fastq.gz

pre-assembled/
└── isolate_contigs.fasta

pacbio/
└── isolate_contigs.fasta
```
**Typical Usage Examples:**
```bash
# Illumina/Illumina+Nanopore data for each isolate
python ardetype.py -i <batch_path> -m all

# Nanopore data
python ardetype.py -i <batch_path> -m all -ont

# Pre-assembled/PacBio (_contigs.fasta) data
python ardetype.py -i <batch_path> -m all -fa

# Merge batches (Illimina + Illumina; Nanopore + Nanopore; Hybrid + Hybrid; Fasta + Fasta) - add -fa or -ont flag for fasta/Nanopore data
python ardetype.py -m merge --merge_from <batch_path1> <batch_path2> <batch_path3> -o <target_batch_path>

# Rerun all rules
python ardetype.py -m all -i <batch_path> -f

# Rerun specific rules (example) - use rule names as defined in snakefiles
python ardetype.py -m all -i <batch_path> -fr resfinder plasmidfinder

```

Notes:
- It is recommended to use the `screen` utility to run the pipeline in the background, especially for long-running jobs or remote sessions.
[- In merge mode, the folder specified by <batch_path> must contain a marker file indicating the sequencing technology used — either ILL_mark, ONT_mark, or HYB_mark.<br>This file is automatically generated by the pipeline during initial processing.]
- The `-fr` option does not validate whether the specified rules exist in the pipeline.<br>Please double-check rule names (`snakefiles/`) and review the output to ensure that the expected steps are reprocessed correctly.
- Refer to `config_files/json/specific_tool_map.json` for details on each tool and the species it is applicable to.
- By default, live pipeline logs are saved in the `Ardetype/.snakemake/logs/` directory.
- The pipeline will iteratively exclude samples that fail any Snakemake rule in any module more than --retry_times (default: 3).
    - All input and output files associated with excluded samples will be moved to the <batch_path>_failed directory.
    - If all samples are excluded, the pipeline will terminate.
    - Detailed error information can be found in the Ardetype/.snakemake/log/ and <batch_path>/logs/ directories.

**Command-line options:**


```text
  -i, --input
      Path to directory that contains files to be analysed.
      Files are NOT parsed from subfolders.
      Check `config_files/json/module_data/patterns/inputs` for expected input file extensions for the `bact_core` module.

  -c, --config
      Path to a custom pipeline configuration file in YAML format.
      If not provided, a default config file (`config_files/yaml/config_modular_local.yaml`) will be used.
      The configuration defines global parameters and tool-specific settings used during execution.

  -o, --output_dir
      Directory where pipeline outputs will be saved.
      By default, it is set to the same value as `--input`, which is recommended for idempotent re-runs.
      In `merge` mode, this must be explicitly specified as the destination for combined results.

  -r, --retry_times
      Number of times Snakemake will attempt to rerun failed jobs.
      Helps recover from intermittent network or I/O errors. Default is 3.

  -n, --num_jobs
      Maximum number of Snakemake jobs to run in parallel (when executing on HPC).
      Default is 12.

  -cpu, --max_cpus  Maximum number of CPUs Snakemake is allowed to use to run the jobs on     
      local machine. Default is 6.

  -hpc, --hpc
      Enable HPC execution mode.
      With appropriate configuration, this allows submitting each pipeline module
      (via `subscripts/ardetype_jobscript.sh`, `--hpc`, and `-sm` flags)
      or per-sample jobs (via `config_files/yaml/cluster.yaml` and `--hpc`) to a cluster scheduler (PBS by default).
      By default, the pipeline runs locally on the machine where the script is invoked.

  -sm, --submit_modules
      Run each pipeline module as a separate HPC job.
      Snakemake sub-jobs are not submitted to the cluster in this mode.
      Requires manual adjustment of `subscripts/ardetype_jobscript.sh` to control behavior.
      Effective only with `--hpc`.

  -f, --force_all
      Force rerun of all Snakemake rules across all pipeline modules,
      regardless of previous results or timestamps.

  -sp, --skip_packing
      Disable output packing.
      By default, output files are organized into per-sample folders.
      Use this flag to keep a flat output directory structure.

  -ont, --nanopore_only
      Run the pipeline on a folder containing exclusively Nanopore reads (with `_ONT.fastq.gz` extension).
      By default, the pipeline expects either Illumina data (`_001.fastq.gz`) or combined Illumina and Nanopore data.

  -mf MERGE_FROM [MERGE_FROM ...], --merge_from MERGE_FROM [MERGE_FROM ...]
      One or more paths to folders containing ARDETYPE output to be merged.
      Allowed combinations: Hybrid + Hybrid, Illumina + Illumina, Nanopore + Nanopore, Fasta + Fasta.

  -fr FORCE_RULES [FORCE_RULES ...], --force_rules FORCE_RULES [FORCE_RULES ...]
      Specify one or more Snakemake rule names (check `snakefiles/`) to force rerun within the pipeline.


required arguments:
  -m, --run_mode
      Selects pipeline execution mode:

        all    - Run the complete pipeline on raw FASTQ files, including all analysis modules:
                   `bact_core`, `bact_shell`, `bact_tip`, `bact_shape`.

        merge  - Combine outputs from previous ARDETYPE runs (see `--merge_from`),
                   then generate unified reports under `--output_dir`.
```

**Report summary**

| File | Description |
|------|-------------|
| aquamis_qc_report.csv | comparing QC values against threshold defined in [AQUAMIS](https://gitlab.com/bfr_bioinformatics/AQUAMIS/-/blob/master/resources/AQUAMIS_thresholds.json?ref_type=heads) pipeline |
| ardetype_report.csv | QC and truncated species-specific typing information |
| agrvate_report.csv | [agr typing](https://github.com/VishnuRaghuram94/AgrVATE) for S.aureus |
| amrfp_mutation_report.csv | [Detecting AMR markers](https://github.com/ncbi/amr) in bacterial NGS data. |
| chewbbaca_qc_report.csv | cgmlst profile quality using [chewbbaca](https://github.com/B-UMMI/chewBBACA) & Ridom schema |
| ectyper_report.csv | [serotyping](https://github.com/phac-nml/ecoli_serotyping) for E.coli |
| emmtyper_report.csv | [emm-type and emm-cluster](https://github.com/MDU-PHL/emmtyper) prediction for S.pyogenes |
| harmonized_resistance_profile.tsv | combined [resfinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) and [RGI/CARD](https://github.com/arpcard/rgi) predictions |
| hicap_report.csv | [serotyping](https://github.com/scwatts/hicap) for H. Influenzae |
| kleborate_report.csv | [serotyping, resistance and virulence](https://github.com/katholt/Kleborate) profiles for K.pneumoniae |
| kraken2contigs_report.csv | all hits with [kraken2](https://github.com/DerrickWood/kraken2)/[minikraken](https://benlangmead.github.io/aws-indexes/k2#:~:text=Older%20%E2%80%9CMinikraken%E2%80%9D%20indexes) applied to contigs |
| kraken2reads_report.csv | all hits with [kraken2](https://github.com/DerrickWood/kraken2)/[minikraken](https://benlangmead.github.io/aws-indexes/k2#:~:text=Older%20%E2%80%9CMinikraken%E2%80%9D%20indexes) applied to contigs |
| legsta_report.csv | [sequence-based typing (SBT)](https://github.com/tseemann/legsta) for L. Pneumophila |
| lissero_report.csv | [serotyping](https://github.com/MDU-PHL/LisSero) for L.monocytogenes |
| lrefinder_report.csv | [predicting linezolid resistance](https://bitbucket.org/genomicepidemiology/lre-finder) in Enterococci faecalis and E. faecium |
| meningotype_report.csv | [serotyping](https://github.com/MDU-PHL/meningotype) for N.meningitidis |
| mobtyper_contig_summary.csv | [Detection and reconstruction](https://github.com/phac-nml/mob-suite) of mobile genetic elements from bacterial NGS data |
| mobtyper_summary.csv | [Detection](https://github.com/phac-nml/mob-suite) and reconstruction of mobile genetic elements from bacterial NGS data |
| plasmidfinder_summary.csv | [Detecting mobile genetic elements](https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/) in bacterial NGS data (contigs) |
| pointfinder_report.csv | [Detecting AMR-associated mutations](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) in bacterial NGS data (contigs) |
| quast_report.csv | [Quality control of genomes](https://github.com/ablab/quast) assembled from short-read NGS data |
| resfinder_mobtyper_mapping.csv | AMR-MGE link prediction from MOBSuite and Resfinder based on contig overlap |
| resfinder_pheno_table_gathered.csv | phenotypic resistance predictions by [resfinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) |
| seqsero2_report.csv | [serotyping](https://github.com/denglab/SeqSero2) for S.enterica |
| seroba_report.csv | [serotyping](https://github.com/sanger-pathogens/seroba) for S.pneumoniae |
| shigatyper_report.csv | [serotyping](https://github.com/CFSAN-Biostatistics/shigatyper) for Shigella |
| sistr_report.csv | [serotyping](https://github.com/phac-nml/sistr_cmd) for S.enterica |
| software_log.csv | software/database versions |
| spatyper_report.csv | [spa typing](https://github.com/HCGB-IGTP/spaTyper) for S.aureus |
| stecfinder_report.csv | [Shigatoxin typing](https://github.com/LanLab/STECFinder) for E.coli |
| virulencefinder_summary.csv | [Detecting virulence genetic elements](https://bitbucket.org/genomicepidemiology/virulencefinder/src/master/) in bacterial NGS data (including full [VFDB](https://www.mgc.ac.cn/VFs/download.htm) database) |
| pf_abaumannii_report.csv | [Acinetobacter baumannii MGA search results](https://github.com/MehradHamidian/AcinetobacterPlasmidTyping) |