# Descriptive and threshold-detection recalculation

The tables use the paired CMAQ–AQS fire-season dataset for June–October 2008–2018. The primary analysis retains 187,051 rows and quarantines 18 rows with CMAQ total PM₂.₅ above 557 µg/m³ without capping. The raw sensitivity retains all 187,069 otherwise-valid rows.

Threshold detection treats an observation as an event when observed PM₂.₅ is at least 12 or 35 µg/m³ and treats CMAQ as detecting the event when CMAQ total PM₂.₅ reaches the same threshold. The 9 µg/m³ annual NAAQS is not used as a daily paired-observation threshold.

## Principal results

- At 12 µg/m³ under primary QC, sensitivity is 0.454, specificity 0.891, precision 0.661, balanced accuracy 0.672, F1 0.539, and frequency bias 0.687. CMAQ therefore misses more than half of observed ≥12 µg/m³ events and underpredicts their frequency.
- At 35 µg/m³, sensitivity is 0.408 and specificity 0.987. Precision is only 0.266, F1 is 0.322, and frequency bias is 1.539. CMAQ produces more modeled ≥35 µg/m³ events than observed events, yet still misses 1,273 of the 2,151 observed exceedances because the modeled and observed events often occur on different paired site-days.
- The primary all-year descriptive result is mean bias −1.40 µg/m³, MAE 5.16 µg/m³, RMSE 9.96 µg/m³, Pearson r 0.502, and r² 0.252.
- In the high observed-concentration regime, primary-QC mean error is −9.79 µg/m³, median error −17.69 µg/m³, MAE 31.94 µg/m³, and RMSE 50.88 µg/m³. Restoring the quarantined extremes changes these to −3.34, −17.57, 38.22, and 150.37 µg/m³, respectively.

## Reconciliation issue

The Rim Fire values reproduce the reported supplementary table after rounding. The current paired file also reproduces the Klamath r² of 0.66, but it yields MAE 12.67, RMSE 23.33, mean bias 8.56 µg/m³, and NMB 45.46%, rather than the reported 15.75, 29.37, 13.04 µg/m³, and 74.0%. No quarantined extreme occurs in 2008, so the discrepancy is unrelated to the new extreme-value rule. The manuscript should use the recalculated values unless an earlier Klamath-specific input dataset or filtering rule can be recovered and justified.

