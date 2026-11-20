# Statistical analysis at 300 Hz

The reanalysis computes Welch power at 4--7 Hz (theta), 12--15 Hz (SMR) and 20--30 Hz (high-beta) using the original 512-sample Hann window and 256-sample overlap. Delta-phase--SMR-amplitude coupling uses Hilbert analytic signals and 200 circular-shift surrogates.

For each EEG feature, the between-group contrast is the trained minus control difference in mean post--pre change. Two-sided group-label permutations enumerate 12,870 balanced assignments conditional on exchangeability. The four primary contrasts are BH-adjusted. Exploratory high-beta sensitivity includes leave-one-trained-participant-out deletion and relative power (20--30 over 0.5--30 Hz), outside the four-primary-outcome FDR family.

The four-feature PCA standardizes trained-group change scores. Bootstraps refit standardization/PCA and align components to the full-sample solution; a control-reference scaling check is also reported. Exact PERMANOVA and dispersion analyses use pooled standardized feature changes. The related manuscript describes study-design limitations.

The original study uses a nonfeedback control, not sham or yoked feedback. This is an observational group-contrast analysis; group-allocation randomization has not been verified.
