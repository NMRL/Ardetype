# ARDETYPE

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/release/NMRL/Ardetype?color=Green)

> **ARDETYPE** is a modular, scalable pipeline for bacterial whole-genome sequencing (WGS) analysis. It analyses **Illumina (paired-end)** reads in `fastq.gz` format and **pre-assembled genomes** in FASTA format (e.g. PacBio assemblies).<br>Designed for clinical microbiology and public health applications, it enables streamlined and reproducible analysis of bacterial isolates through both species-agnostic and species-specific workflows.


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
| [AcinetobacterPlasmidTyping](https://github.com/MehradHamidian/AcinetobacterPlasmidTyping) | bact_tip   |
| [agrvate](https://github.com/VishnuRaghuram94/AgrVATE)                                     | bact_tip   |
| [amrfinder+](https://github.com/ncbi/amr)                                                  | bact_shell |
| [blast+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)             | bact_tip   |
| [capybara](https://github.com/Zhou-lab-SUDA/CAPYBARA)                                      | bact_tip   |
| [cgmlst-dists](https://github.com/tseemann/cgmlst-dists)                                   | bact_shape |
| [chewbbaca](https://github.com/B-UMMI/chewBBACA)                                           | bact_tip   |
| [ectyper](https://github.com/phac-nml/ecoli_serotyping)                                    | bact_tip   |
| [el_gato](https://github.com/CDCgov/el_gato)                                               | bact_tip   |
| [emmtyper](https://github.com/MDU-PHL/emmtyper)                                            | bact_tip   |
| [fastp](https://github.com/OpenGene/fastp)                                                 | bact_core  |
| [hamronization](https://github.com/pha4ge/hAMRonization)                                   | bact_shell |
| [hicap](https://github.com/scwatts/hicap)                                                  | bact_tip   |
| [kleborate](https://github.com/katholt/Kleborate)                                          | bact_tip   |
| [krakentools](https://github.com/jenniferlu717/KrakenTools)                                | bact_core  |
| [kraken2](https://github.com/DerrickWood/kraken2)                                          | bact_core  |
| [krona](https://github.com/marbl/Krona)                                                    | bact_core  |
| [lissero](https://github.com/MDU-PHL/LisSero)                                              | bact_tip   |
| [lrefinder](https://bitbucket.org/genomicepidemiology/lre-finder)                          | bact_tip   |
| [meningotype](https://github.com/MDU-PHL/meningotype)                                      | bact_tip   |
| [mlst](https://github.com/tseemann/mlst)                                                   | bact_shell |
| [mob-suite](https://github.com/phac-nml/mob-suite)                                         | bact_shell |
| [ngmaster](https://github.com/MDU-PHL/ngmaster)                                            | bact_tip   |
| [pigz](https://zlib.net/pigz/)                                                             | bact_core  |
| [plasmidfinder](https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/)       | bact_shell |
| [quast](https://github.com/ablab/quast)                                                    | bact_shell |
| [resfinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/)               | bact_shell |
| [rgi](https://github.com/arpcard/rgi)                                                      | bact_shell |
| [seqsero2](https://github.com/denglab/SeqSero2)                                            | bact_tip   |
| [seroba](https://github.com/sanger-pathogens/seroba)                                       | bact_tip   |
| [shigatyper](https://github.com/CFSAN-Biostatistics/shigatyper)                            | bact_tip   |
| [shovill](https://github.com/tseemann/shovill)                                             | bact_core  |
| [sistr](https://github.com/phac-nml/sistr_cmd)                                             | bact_tip   |
| [spatyper](https://github.com/HCGB-IGTP/spaTyper)                                          | bact_tip   |
| [staphopia-sccmec](https://github.com/staphopia/staphopia-sccmec)                          | bact_tip   |
| [stecfinder](https://github.com/LanLab/STECFinder)                                         | bact_tip   |
| [virulencefinder](https://bitbucket.org/genomicepidemiology/virulencefinder/src/master/)   | bact_shell |
| [rmlst API](https://pubmlst.org/species-id/species-identification-via-api)                 | bact_shell |
| [PLSDB API](https://github.com/CCB-SB/plsdbapi)                                            | bact_shell |
---

## Installation

```bash
git clone -b master https://github.com/NMRL/Ardetype
cd Ardetype
conda env create -f ardetype.yaml
```

- Note: Currently PLSDB API module requires [some adjustments](https://github.com/CCB-SB/plsdbapi/issues/4) to the source code in order to make it work.
- [Finding python module files in conda environment](https://stackoverflow.com/questions/31003994/where-is-site-packages-located-in-a-conda-environment).

### Container images and databases

The Singularity images and databases are published as two archives in one Zenodo record: **Ardetype: Singularity images and reference databases**, https://doi.org/10.5281/zenodo.24167460

| Archive | Contents | Size |
|---|---|---|
| `ardetype_images_2026-11-12.tar.gz` | 35 Singularity images: 29 rebuilt from conda-forge/bioconda lock files and 6 original builds (PlasmidFinder, MOB-suite, AgrVATE, staphopia-sccmec, ngmaster, BLAST+), plus the recipe to rebuild the 29 | ≈18 GB |
| `ardetype_databases_2026-11-12.tar.gz` | Kraken2 (PhiX, human T2T, MiniKraken), ResFinder, VirulenceFinder + VFDB, SeroBA, LRE-Finder, PlasmidFinder, *A. baumannii* plasmid typing, AQUAMIS QC thresholds, 21 cgMLST schemes | ≈12 GB (≈42 GB extracted) |

All paths in the archives start at the repository root, so each archive is extracted there in one step. Besides the images or databases, each archive contains `resources/release_docs/<images|databases>_2026-11-12/` with a `README.md`, `references.txt` (Harvard-style citations of every tool and database; please cite them), `MANIFEST.tsv` and `checksums/<item>.sha256` (SHA-256 of every file, one checksum file per image or database).

**1. Download** both archives from the [Zenodo record](https://doi.org/10.5281/zenodo.24167460) and compare each with the MD5 checksum Zenodo lists for it:
```bash
md5sum ardetype_images_2026-11-12.tar.gz ardetype_databases_2026-11-12.tar.gz
```

**2. Extract** in the repository root:
```bash
tar -xzf <download folder>/ardetype_images_2026-11-12.tar.gz    -C .
tar -xzf <download folder>/ardetype_databases_2026-11-12.tar.gz -C .
```
Use plain `tar -xzf`: the archives keep file modification times, which the pipeline uses to report database dates.

**3. Verify the extracted files** (in the repository root):
```bash
sha256sum --quiet -c resources/release_docs/images_2026-11-12/checksums/*.sha256    && echo "images OK"
sha256sum --quiet -c resources/release_docs/databases_2026-11-12/checksums/*.sha256 && echo "databases OK"
```

**Scope of the archives:** everything needed to run all four modules (`bact_core`, `bact_shell`, `bact_tip`, `bact_shape`) on Illumina paired-end reads and on assemblies (`-fa`). The *Mycobacterium* cgMLST scheme is not included and its entry in `config_files/yaml/config_modular_local.yaml` is disabled.

The databases built into images (AMRFinderPlus database, CARD, PubMLST MLST schemes, MOB-suite plasmid databases) come with the images archive. rMLST and PLSDB are online services (see Configuration).

To rebuild the 29 conda-forge/bioconda images instead of downloading them, see `resources/release_docs/images_2026-11-12/build_recipe/README.md`.

---

## Configuration

### Singularity

All tools run in Singularity images. Set the Singularity executable explicitly in `config_files/yaml/config_modular_local.yaml`:

```yaml
singularity_bin:
  "/usr/bin/singularity"
```

- Find the path with `command -v singularity` (or `command -v apptainer`) on a node where Singularity is available, **outside** the activated `ardetype` conda environment. If your system provides it through environment modules, run `module load <singularity module>` first.
- The path must be the same on every node that runs pipeline jobs (login node and compute nodes in HPC mode).
- The pipeline links this executable as `singularity` at the front of `PATH`, so Snakemake and every job it submits use it. It prints `Using singularity: <path>` at start.
- If `singularity_bin` is empty, the first `singularity` on `PATH` is used. In the activated conda environment this is the `singularityce` package from `ardetype.yaml`. It works, but it cannot use the setuid starter and unpacks every image into `SINGULARITY_TMPDIR` before each run, which is much slower and needs a lot of temporary space. The pipeline prints a warning in that case.
- Image unpacking and cache directories default to `tmp/` in the repository (`SINGULARITY_TMPDIR`, `SINGULARITY_CACHEDIR` and their `APPTAINER_*` equivalents); set these variables before starting the pipeline to use another location.

### rMLST

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

### Other settings

The `config_files/yaml/config_modular_local.yaml` file can be customized—for example, to define a custom location for the scratch directory used for intensive I/O operations.

- Relative paths should be specified with respect to the root of the `Ardetype` repository folder.
- Absolute paths are also supported, as long as they are accessible from the machine running the pipeline.

---

## Usage

**Input Structure Examples:**

- `isolate` stands for sample ID.
- Having `_` as part of sample ID may cause unexpected behavior.
```bash
illumina/
├── isolate_R1_001.fastq.gz
└── isolate_R2_001.fastq.gz

pre-assembled/          # any assembly, e.g. from PacBio reads
└── isolate_contigs.fasta
```
**Typical Usage Examples:**
```bash
# Illumina paired-end data
python ardetype.py -i <batch_path> -m all

# Pre-assembled (_contigs.fasta) data
python ardetype.py -i <batch_path> -m all -fa

# Merge batches (Illumina + Illumina; Fasta + Fasta) - add the -fa flag for fasta data
python ardetype.py -m merge --merge_from <batch_path1> <batch_path2> <batch_path3> -o <target_batch_path>

# Rerun all rules
python ardetype.py -m all -i <batch_path> -f

# Rerun specific rules (example) - use rule names as defined in snakefiles
python ardetype.py -m all -i <batch_path> -fr resfinder plasmidfinder

```

Notes:
- It is recommended to use the `screen` utility to run the pipeline in the background, especially for long-running jobs or remote sessions.
- In merge mode, the folder specified by <batch_path> must contain a marker file indicating the sequencing technology used — either ILL_mark (Illumina) or FA_mark (assemblies).<br>This file is automatically generated by the pipeline during initial processing.
- The `-fr` option does not validate whether the specified rules exist in the pipeline.<br>Please double-check rule names (`snakefiles/`) and review the output to ensure that the expected steps are reprocessed correctly.
- Refer to `config_files/json/specific_tool_map.json` for details on each tool and the species it is applicable to.
- By default, live pipeline logs are saved in the `Ardetype/.snakemake/log/` directory.
- The pipeline will iteratively exclude samples that fail any Snakemake rule in any module more than --retry_times (default: 3).
    - All input and output files associated with excluded samples will be moved to the <batch_path>_failed directory.
    - If all samples are excluded, the pipeline will terminate.
    - Detailed error information can be found in the Ardetype/.snakemake/log/ and <batch_path>/logs/ directories.

**Command-line options:**


```text
  -i, --input
      Path to directory that contains files to be analysed.
      Files are NOT parsed from subfolders.
      Expected file name endings: Illumina `_R1_001.fastq.gz` / `_R2_001.fastq.gz`,
      assemblies `_contigs.fasta`
      (defined in `config_files/json/module_data.json`, section "core" -> "patterns" -> "inputs").

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
      or per-sample jobs (via `config_files/yaml/cluster_slurm.yaml` and `--hpc`) to a cluster scheduler.
      The scheduler is set by the `cluster:` key of that file (`slurm` or `pbs`).
      By default, the pipeline runs locally on the machine where the script is invoked.

  -sm, --submit_modules
      Run each pipeline module as a separate HPC job.
      Snakemake sub-jobs are not submitted to the cluster in this mode.
      Works on PBS/Torque only (modules are submitted with `qsub -F` and polled with `qstat`).
      Not usable on Slurm: there, use `--hpc` without `-sm` (per-rule jobs via `cluster_slurm.yaml`).
      Requires manual adjustment of `subscripts/ardetype_jobscript.sh` to control behavior.
      Effective only with `--hpc`.

  -f, --force_all
      Force rerun of all Snakemake rules across all pipeline modules,
      regardless of previous results or timestamps.

  -sp, --skip_packing
      Disable output packing.
      By default, output files are organized into per-sample folders.
      Use this flag to keep a flat output directory structure.

  -fa, --fasta
      Run the pipeline on a folder containing exclusively assemblies, one per isolate,
      named `<isolate>_contigs.fasta` (e.g. pre-assembled genomes or PacBio assemblies).
      Read-based steps (QC, read filtering, assembly) are skipped.
      By default, the pipeline expects Illumina paired-end data (`_R1_001.fastq.gz` / `_R2_001.fastq.gz`).

  -mf MERGE_FROM [MERGE_FROM ...], --merge_from MERGE_FROM [MERGE_FROM ...]
      One or more paths to folders containing ARDETYPE output to be merged.
      Allowed combinations: Illumina + Illumina, Fasta + Fasta.

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
| YYYY-MM-DD_aquamis_qc_report.csv | comparing QC values against threshold defined in [AQUAMIS](https://gitlab.com/bfr_bioinformatics/AQUAMIS/-/blob/master/resources/AQUAMIS_thresholds.json?ref_type=heads) pipeline |
| YYYY-MM-DD_ardetype_report.csv | QC and truncated species-specific typing information (prefixed with the date of the run) |
| agrvate_report.csv | [agr typing](https://github.com/VishnuRaghuram94/AgrVATE) for S.aureus |
| amrfp_mutation_report.csv | status of every resistance-associated point mutation position screened by [AMRFinderPlus](https://github.com/ncbi/amr) (`--mutation_all`) in the assembled contigs, for supported organisms. **Warning:** AMRFinderPlus accepts only assemblies, so each position reports the consensus base of the assembly; minority variants (e.g. heteroresistance, mixed populations or variants in some copies of multi-copy genes such as 23S rRNA) are not detected |
| amrfp_report.csv | AMR, stress response and virulence genes and point mutations found by [AMRFinderPlus](https://github.com/ncbi/amr) |
| capybara_report.csv | [ESL lineage and clade](https://github.com/Zhou-lab-SUDA/CAPYBARA) detection for A. baumannii |
| chewbbaca_qc_report.csv | cgMLST profile quality using [chewbbaca](https://github.com/B-UMMI/chewBBACA) with schemes from [cgMLST.org](https://www.cgmlst.org/ncs) (Ridom), [EnteroBase](https://enterobase.warwick.ac.uk/), [PubMLST](https://pubmlst.org/) and [Zenodo](https://zenodo.org/records/17871974) (LIT *L. pneumophila* cg/wgMLST) |
| cgmlst_aggregated/ | cgMLST allele profiles from [chewbbaca](https://github.com/B-UMMI/chewBBACA), one `agg_results_alleles_<species>.tsv` per species |
| ectyper_report.csv | [serotyping](https://github.com/phac-nml/ecoli_serotyping) for E.coli |
| emmtyper_report.csv | [emm-type and emm-cluster](https://github.com/MDU-PHL/emmtyper) prediction for S.pyogenes |
| harmonized_resistance_profile.tsv | combined [resfinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) and [RGI/CARD](https://github.com/arpcard/rgi) predictions |
| hicap_report.csv | [serotyping](https://github.com/scwatts/hicap) for H. Influenzae |
| kleborate_report.csv | [serotyping, resistance and virulence](https://github.com/katholt/Kleborate) profiles for K.pneumoniae |
| kraken2contigs_report.csv | all hits with [kraken2](https://github.com/DerrickWood/kraken2)/[minikraken](https://benlangmead.github.io/aws-indexes/k2#:~:text=Older%20%E2%80%9CMinikraken%E2%80%9D%20indexes) applied to contigs |
| kraken2reads_report.csv | all hits with [kraken2](https://github.com/DerrickWood/kraken2)/[minikraken](https://benlangmead.github.io/aws-indexes/k2#:~:text=Older%20%E2%80%9CMinikraken%E2%80%9D%20indexes) applied to reads |
| elgato_report.csv | [sequence-based typing (SBT)](https://github.com/CDCgov/el_gato) for L. Pneumophila |
| lissero_report.csv | [serotyping](https://github.com/MDU-PHL/LisSero) for L.monocytogenes |
| lrefinder_report.csv | [predicting linezolid resistance](https://bitbucket.org/genomicepidemiology/lre-finder) in Enterococci faecalis and E. faecium |
| mlst_long_report.csv | [MLST](https://github.com/tseemann/mlst) alleles, one row per locus |
| mlst_raw_report.csv | [MLST](https://github.com/tseemann/mlst) scheme, sequence type (ST) and allelic profile, one row per sample |
| meningotype_report.csv | [serotyping](https://github.com/MDU-PHL/meningotype) for N.meningitidis |
| mobtyper_contig_summary.csv | [Detection and reconstruction](https://github.com/phac-nml/mob-suite) of mobile genetic elements from bacterial NGS data |
| mobtyper_summary.csv | [Detection](https://github.com/phac-nml/mob-suite) and reconstruction of mobile genetic elements from bacterial NGS data |
| plasmidfinder_summary.csv | [Detecting mobile genetic elements](https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/) in bacterial NGS data (contigs) |
| pointfinder_report.csv | [Detecting AMR-associated mutations](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) in bacterial NGS data (contigs) |
| quast_report.csv | [Quality control of genomes](https://github.com/ablab/quast) assembled from short-read NGS data |
| resfinder_mobtyper_mapping.csv | AMR-MGE link prediction from MOBSuite and Resfinder based on contig overlap |
| resfinder_reads_report.csv | acquired AMR genes found by [resfinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) in the reads |
| rmlst_loci_report.csv | species, support and ribosomal ST from the [rMLST](https://pubmlst.org/species-id) API, with rMLST alleles |
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
| plsdb_report.csv | [PLSDB similarity search results for MOB-Suite constructed plasmid sequences](https://github.com/CCB-SB/plsdbapi) |