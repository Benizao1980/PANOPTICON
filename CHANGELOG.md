# Changelog

## 0.2.0 - proposed

### Added
- tool-neutral `pangenome_ecology.py`
- observed and rarefied group-sharing metrics
- group-core overlap and group-specific family summaries
- exact nucleotide/protein allele-sharing input
- robust Panaroo Rtab-only loading
- metadata preparation and host/ecology grouping
- complete animal-associated *S. aureus* tutorial
- example Slurm jobs and metadata template
- smoke tests and GitHub Actions workflow

### Changed
- renamed generic host stratification workflow to `pangenome_stratification.py`
- clarified that rarefied metrics should be preferred over raw overlap
- documented memory requirements for large datasets

### Compatibility
- legacy `--min-host-count` and `--min-host-prevalence` aliases remain accepted by `pangenome_ecology.py`
