# ARDETYPE

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/release/NMRL/Ardetype?color=Green)

**ARDETYPE** is a modular pipeline for processing whole genome sequencing (WGS) data from both **Illumina (paired-end)** and **Oxford Nanopore** reads. It supports both species-agnostic and species-specific analyses.

The pipeline is structured into four Snakemake modules, each representing a distinct phase of the analysis:

- `bact_core`: Core tasks such as quality control (QC), contamination filtering, and species identification.
- `bact_shell`: General-purpose tools for resistance gene detection, plasmid analysis, and more.
- `bact_tip`: Species-specific typing tools.
- `shape`: Aggregation and result integration.

## Tools by Module

| Tool                                                                 | Module     |
|----------------------------------------------------------------------|------------|
| [agrvate](https://github.com/VishnuRaghuram94/AgrVATE)              | bact_tip   |
| [amrfinder+](https://github.com/ncbi/amr)                            | bact_shell |
| [capybara](https://github.com/Zhou-lab-SUDA/CAPYBARA)                | bact_tip   |
| [chewbbaca](https://github.com/B-UMMI/chewBBACA)                     | bact_tip   |
| [circlator](https://github.com/sanger-pathogens/circlator)          | bact_core  |
| [ectyper](https://github.com/phac-nml/ecoli_serotyping)             | bact_tip   |
| [emmtyper](https://github.com/MDU-PHL/emmtyper)                      | bact_tip   |
| [fastp](https://github.com/OpenGene/fastp)                           | bact_core  |
| [filtlong](https://github.com/rrwick/Filtlong)                       | bact_core  |
| [flye](https://github.com/fenderglass/Flye)                          | bact_core  |
| [hicap](https://github.com/scwatts/hicap)                            | bact_tip   |
| [kleborate](https://github.com/katholt/Kleborate)                    | bact_tip   |
| [kraken2](https://github.com/DerrickWood/kraken2)                    | bact_core  |
| [legsta](https://github.com/tseemann/legsta)                         | bact_tip   |
| [lissero](https://github.com/MDU-PHL/LisSero)                        | bact_tip   |
| [lrefinder](https://bitbucket.org/genomicepidemiology/lre-finder)   | bact_tip   |
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
| [shigatyper](https://github.com/CFSAN-Biostatistics/shigatyper)     | bact_tip   |
| [shovill](https://github.com/tseemann/shovill)                       | bact_core  |
| [sistr](https://github.com/phac-nml/sistr_cmd)                       | bact_tip   |
| [snikt](https://github.com/piyuranjan/SNIKT)                         | bact_core  |
| [spatyper](https://github.com/HCGB-IGTP/spaTyper)                    | bact_tip   |
| [staphopia-sccmec](https://github.com/staphopia/staphopia-sccmec)   | bact_tip   |
| [stecfinder](https://github.com/LanLab/STECFinder)                   | bact_tip   |
| [virulencefinder](https://bitbucket.org/genomicepidemiology/virulencefinder/src/master/) | bact_shell |

---

## Installation

```bash
git clone -b standalone_rework https://github.com/NMRL/Ardetype
cd Ardetype
conda env create -f ardetype.yaml

# Downloading pre-packaged resources (databases & singularity image files)
rm -rf resources
wget <archive_url>
tar -xvzf resources.tar.gz
```
## Configuration

To enable **[rMLST](https://pubmlst.org/species-id)** support, you need to configure API access to the PubMLST REST service. This requires creating a `config_files/.env` file with authentication credentials.

```env
CONSUMER_KEY=""
CONSUMER_SECRET=""
TEST_REST_URL="https://rest.pubmlst.org/db/pubmlst_rmlst_seqdef_kiosk/schemes/1"
SESSION_REST_URL="https://rest.pubmlst.org/db/pubmlst_test_seqdef"
TEST_WEB_URL="https://pubmlst.org/bigsdb?db=pubmlst_test_seqdef"
ACCESS_TOKEN_PATH="<path_for_access_token>"
```

- You must request an API key pair from PubMLST to populate the `CONSUMER_KEY` and `CONSUMER_SECRET`.
- See the [official documentation](https://bigsdb.readthedocs.io/en/latest/rest.html#authentication) for details on how to do this.
- Replace `CONSUMER_KEY` and `CONSUMER_SECRET` with your credentials from PubMLST.
- `ACCESS_TOKEN_PATH` should point to a local file where the access token will be saved automatically after authentication.
- In case rmlst configuration was not done, species will be inferred based on kraken2/minikraken contig classfication.
- **Note**: If rMLST-based species identification is not configured, the species will be inferred from contig classification using Kraken2 with the MiniKraken database.

## Usage

- Input Requirements:

     - Illumina: Paired-end FASTQ files named as <sample_id>_R1_001.fastq.gz and <sample_id>_R2_001.fastq.gz
     - Nanopore: Single-end FASTQ file named as <sample_id>_ONT.fastq.gz

- Run for Illumina or Illumina + Nanopore:

     ```python ardetype.py -i <batch_folder> -o <batch_folder> -m all -rl -p```

- Run for Nanopore only:

     ```python ardetype.py -i <batch_folder> -o <batch_folder> -m all -rl -p -ont```

