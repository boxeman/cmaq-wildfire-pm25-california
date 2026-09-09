# CMAQ wildfire PM2.5 evaluation in California

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22677978.svg)](https://doi.org/10.5281/zenodo.22677978)

This repository contains the frozen data, code, fold assignments, diagnostic outputs, and recalculated tables supporting the manuscript *Evaluating CMAQ Representation of Wildfire PM2.5 in California (2008–2018): Regime-Dependent Performance and Nonlinear Error Structure*.

**Archive status:** GitHub release `v1.0.0` is archived in Zenodo at [https://doi.org/10.5281/zenodo.22677978](https://doi.org/10.5281/zenodo.22677978). The all-versions concept DOI is [https://doi.org/10.5281/zenodo.22677977](https://doi.org/10.5281/zenodo.22677977).

## Study scope

The analysis evaluates daily paired Community Multiscale Air Quality (CMAQ) estimates and monitoring observations during California fire seasons (June–October), 2008–2018. It includes conventional descriptive metrics, paired event detection at 12 and 35 µg m⁻³, Klamath and Rim Fire case studies, and a grouped leakage-resistant Random Forest diagnostic of CMAQ error.

The primary quality-control population contains 187,051 site-days. Eighteen otherwise-valid rows with CMAQ total PM2.5 above 557 µg m⁻³ are quarantined without capping and retained only in explicitly labelled raw-extreme sensitivity analyses.

## Repository contents

- `data/paired_cmaq_aqs_fireseason_2008_2018.csv`: frozen paired input dataset.
- `leakage_resistant_rerun/`: notebook, Python implementation, pinned environment, deterministic site-year assignments, untouched-test predictions, model-comparison tables, interpretation data, and replacement figures.
- `descriptive_threshold_tables/`: recalculation script and frozen descriptive, regime, case-study, threshold-detection, QC, and reconciliation outputs.
- `CITATION.cff`: machine-readable citation metadata for GitHub and Zenodo.
- `MANIFEST.sha256`: SHA-256 checksums for the archived package.

## Reproduction

Python 3.12 was used for the frozen run. Create a clean environment and install the pinned packages:

```bash
python -m venv .venv
python -m pip install -r leakage_resistant_rerun/requirements.txt
```

To rerun the leakage-resistant analysis, open `leakage_resistant_rerun/Bonilla_leakage_resistant_CMAQ_rerun.ipynb`. The notebook has `RUN_FULL_ANALYSIS = True`, resolves the repository data file automatically, and writes regenerated files to `leakage_resistant_rerun/reproduction_check/`.

To recalculate the descriptive and threshold-detection CSV tables:

```bash
python descriptive_threshold_tables/recalculate_descriptive_threshold_tables.py
```

The frozen input SHA-256 is `f2e278ad52373f3a1a6f3d5b01f20f38d9c948d261ea34fe7e4941a33265070d`. The executed notebook and independent comparison reproduced the deterministic fold assignments, ten CSV outputs, and six figures; details are recorded in `leakage_resistant_rerun/reproduction_verification.json`.

## Principal frozen results

- Primary paired analysis: mean bias −1.40 µg m⁻³, MAE 5.16 µg m⁻³, RMSE 9.96 µg m⁻³, Pearson *r* 0.502, and *r*² 0.252.
- Event detection at 12 µg m⁻³: sensitivity 0.454, specificity 0.891, precision 0.661, and frequency bias 0.687.
- Event detection at 35 µg m⁻³: sensitivity 0.408, specificity 0.987, precision 0.266, and frequency bias 1.539.
- Leakage-resistant untouched test: *R*² 0.711, RMSE 5.82 µg m⁻³, and MAE 3.44 µg m⁻³ across 37,430 observations in 428 site-year groups.
- Restoring the 18 quarantined records reduced untouched-test *R*² to 0.381 and increased RMSE to 9.44 µg m⁻³.

## Citation

Please cite the archived version used for the manuscript:

> Perez-Chavez, J., et al. (2026). *Reproducibility package for evaluating CMAQ representation of wildfire PM2.5 in California 2008–2018* (Version 1.0.0). Zenodo. https://doi.org/10.5281/zenodo.22677978

The version DOI identifies the exact frozen release. The concept DOI, `10.5281/zenodo.22677977`, resolves to the latest archived version.

## Data provenance

The paired file contains daily monitoring observations and CMAQ-derived fields used in the manuscript analysis. Monitoring data originated from U.S. EPA air-quality monitoring networks, and the paired CMAQ–AQS dataset was processed with the Atmospheric Model Evaluation Tool. Users should also cite the original data providers and the associated manuscript.
