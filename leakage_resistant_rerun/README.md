# Leakage-resistant CMAQ rerun package

This package independently reconstructs the analysis requested after the leakage concern. It uses the paired CMAQ–AQS fire-season dataset in the manuscript archive and does not depend on a deliverable from José.

## Final analysis design

- Target: `cmaq_total_pm25 - observed_pm25`.
- Group unit: monitoring site × year.
- Test set: deterministic 20% holdout of complete site-year groups (428 groups; 37,430 primary-analysis rows).
- Model selection assessment: five-fold grouped cross-validation within training groups only.
- Primary predictors: CMAQ fire-only PM₂.₅, latitude, longitude, year, month, day of year, distance to fire plus its missingness indicator, and daily observation count.
- Excluded from the primary feature matrix: observed PM₂.₅, total CMAQ PM₂.₅, no-fire CMAQ PM₂.₅, observed regime labels, fire names, and fire-window tags.
- QC: 18 rows with total CMAQ PM₂.₅ > 557 µg/m³ are quarantined from the primary analysis, never capped, and retained in the raw-data sensitivity analysis.

## Headline results

The primary leakage-resistant random forest obtained R² = 0.711, RMSE = 5.82 µg/m³, and MAE = 3.44 µg/m³ on the untouched grouped test set. The linear baseline obtained R² = 0.600 and the fire-only random forest R² = 0.592. Five grouped CV folds ranged from R² = 0.453 to 0.647 (mean 0.581), showing meaningful spatial-temporal heterogeneity.

When the 18 quarantined extremes were restored without capping, test R² fell to 0.381 and RMSE rose to 9.44 µg/m³. This should be reported as a raw-data sensitivity result, not pooled with the primary estimate.

Observed-regime stratification is post hoc. In the high observed-PM₂.₅ regime (≥35 µg/m³; n=473 test rows), RMSE was 29.23 µg/m³ and the model understated the magnitude of CMAQ underprediction by 16.97 µg/m³ on average. The manuscript should therefore avoid presenting the overall R² as evidence of reliable correction during the most extreme smoke conditions.

## Contents

- `Bonilla_leakage_resistant_CMAQ_rerun.ipynb`: readable, rerunnable notebook entry point.
- `leakage_resistant_pipeline.py`: complete implementation.
- `fold_assignments.csv`: source-row split, CV fold, group, and QC status.
- `grouped_cv_metrics.csv`, `grouped_cv_summary.csv`, and `grouped_cv_oof_predictions.csv`: grouped CV evidence.
- `model_comparison_test_metrics.csv`: untouched-test primary, baseline, and sensitivity comparisons.
- `primary_test_metrics_by_regime.csv`: post-hoc regime metrics.
- `Table_S3_regime_error_summary.csv`: primary-QC and raw-sensitivity descriptive summaries.
- `untouched_test_predictions.csv`: row-level held-out predictions.
- `permutation_importance_test_set.csv` and `ale_values.csv`: replacement interpretation data.
- `run_manifest.json`: feature exclusions, QC rules, seed, parameters, package versions, and source SHA-256.
- Replacement figures: Figures 6a, 6b, S10, S11, S12, and S13.
- `manuscript_replacement_text.md` and `replacement_figure_captions.md`: insertion-ready draft language.


## Reproduction

Reproduction was completed in a clean virtual environment created from `requirements.txt`. The notebook has `RUN_FULL_ANALYSIS = True`, the default `DATA_PATH` resolves to the repository's `data` directory, and regenerated outputs are written to `reproduction_check` rather than over the frozen canonical files. The executed notebook's internal comparison passed: fold assignments were identical and seven core metric/interpretation tables matched to a numeric tolerance of 1×10⁻¹². Independent verification extended this to 10 CSV outputs, and all six replacement figures were exact SHA-256 matches. See `reproduction_verification.json` for the recorded result.
