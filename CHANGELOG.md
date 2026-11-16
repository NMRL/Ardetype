# Changelog

## v1.0.1 (2026-11-16)

Container images and databases: archives dated 2026-11-12 in the Zenodo record
[Ardetype: Singularity images and reference databases](https://doi.org/10.5281/zenodo.24167460)
(`ardetype_images_2026-11-12.tar.gz`, `ardetype_databases_2026-11-12.tar.gz`).

### Images and databases
- 29 images rebuilt from conda-forge/bioconda lock files (micromamba base); 6 original builds kept as is.
  Results can differ slightly from earlier image builds.
- `pigz` image updated to 2.8.
- Images and databases distributed as two archives with per-item SHA-256 checksums, `MANIFEST.tsv` and
  `references.txt`; the recipe to rebuild the 29 images is included in the images archive.
- *Mycobacterium* cgMLST scheme not included; its config entry is disabled. *Clostridium* scheme re-enabled.
- cgMLST locus counts corrected (*S. enterica* 2999, *P. aeruginosa* 3866).

### Configuration
- New `singularity_bin` setting: path of the Singularity/Apptainer executable to use (see README, Configuration).
  If empty, the first `singularity` on `PATH` is used, with a warning when this is the conda build.

### Fixes
- Local mode: Snakemake failures were ignored (`snakemake | tee` without `pipefail`), so a failed module did not stop
  the run and failed samples were not moved to `<batch>_failed`. Local mode now also recognises
  "output files are missing" failures.
- PointFinder did not run in Illumina mode since 2025-05-22 (species map no longer lower-cased) and never ran for
  assembly input (species map built after replacing spaces with underscores).
- ResFinder on reads used R1 twice; R2 is now analysed.
- hicap: capsule-positive samples crash while hicap writes its GenBank file (after writing the result table), which
  failed the job and could leave the sample reported as NTHi. The result table is now kept whenever it exists; NTHi is
  reported only when hicap finds no capsule locus genes; other failures are reported as `hicap_failed`.
- SeroBA: "; nan" no longer appended to the serotype.
- Typing columns, `MLST7_ST` and `rMLST` are read as text when building `ardetype_report`
  (numeric-looking types such as SeroBA `3` were shown as `3.0`).
- Removed a leftover `module load singularity` from the assembly-input shell module.

### Documentation
- README: installation from the Zenodo archives, Singularity setup, report descriptions updated
  (including a warning that `amrfp_mutation_report.csv` is based on assemblies and does not detect minority variants),
  `-sm` documented as PBS/Torque only.
