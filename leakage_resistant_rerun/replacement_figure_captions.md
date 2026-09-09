# Replacement figure captions

**Figure 6a.** Observed versus predicted CMAQ total PM₂.₅ error for the untouched leakage-resistant test partition. Error was defined as CMAQ total minus observed PM₂.₅. The test partition comprised 37,430 observations from 428 site-year groups that were not used for model fitting or cross-validation. Hexagon color denotes the logarithm of observation count; the dashed line is the 1:1 relation. Axes are limited to the central 99% of plotted values for readability.

**Figure 6b.** Observed versus predicted CMAQ total PM₂.₅ error in the untouched test partition, stratified post hoc by observed PM₂.₅ concentration regime. Regime labels were not supplied to the model. Panels show background (<12 µg m⁻³), moderate (12–35 µg m⁻³), and high (>35 µg m⁻³) observations. Hexagon color denotes the logarithm of observation count; dashed lines are the 1:1 relation.

**Figure S10.** Permutation importance for the leakage-resistant Random Forest, calculated on a fixed random sample of 20,000 observations from the untouched site-year test partition with 10 repeated permutations. Importance is the decrease in test-set coefficient of determination (R²); error bars denote the standard deviation across repeats. Observed PM₂.₅, total CMAQ PM₂.₅, no-fire CMAQ PM₂.₅, and observed-regime labels were excluded from the primary feature matrix.

**Figure S11.** One-dimensional accumulated local effect (ALE) of CMAQ fire-only PM₂.₅ on predicted CMAQ total PM₂.₅ error in the leakage-resistant model. ALE was estimated from 20 equal-frequency bins in a fixed 20,000-row sample of the untouched test partition and centered at zero. A symmetric logarithmic x-axis is used to display the strongly right-skewed fire-only distribution; the uppermost bin spans the sparse extreme tail and should be interpreted cautiously.

**Figure S12.** One-dimensional accumulated local effect (ALE) of longitude, the next most influential contextual predictor after fire-only PM₂.₅ in test-set permutation analysis. ALE was estimated from 20 equal-frequency bins in a fixed 20,000-row sample of the untouched test partition and centered at zero. This diagnostic describes model behavior and should not be interpreted causally.

**Figure S13.** Distribution of CMAQ total PM₂.₅ error (CMAQ minus observed) by observed concentration regime after application of the primary QC rule. Boxes show the interquartile range and median, whiskers extend to 1.5 times the interquartile range, and diamonds show means. Observed-regime labels were used only for descriptive evaluation, not model training.

