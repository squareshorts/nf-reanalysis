# v2.0.0 — 300-Hz C3 EEG reanalysis and updated reproducibility results

This release provides the analysis implementation and aggregate numerical results supporting the revised manuscript *Responder heterogeneity and multivariate organization of resting-state EEG following SMR neurofeedback* (NeuroImage: Reports).

The 300-Hz source-based analysis supersedes the derived-feature analysis distributed in v1.0.0. The earlier release remains archived as a historical version.

## Analysis included

- Source-based recomputation of 32 pre/post C3 resting-state EEG observations at the 300-Hz acquisition rate.
- Theta (4–7 Hz), SMR (12–15 Hz), and high-beta (20–30 Hz) Welch band powers and surrogate-normalized delta-phase–SMR-amplitude coupling.
- Trained-minus-control changes with exact two-sided permutation inference, bootstrap intervals, effect sizes, and Benjamini–Hochberg adjustment across four EEG outcomes.
- Control-referenced responder classifications and threshold-sensitivity analysis.
- Four-feature PCA, resampling and participant-deletion sensitivity, and exact PERMANOVA/dispersion comparisons.
- Absolute and broadband-normalized high-beta sensitivity analyses.

## Principal aggregate findings

- Absolute high-beta contrast: −1.452, exact p = 0.03403; four-feature FDR q = 0.13613. No primary feature remains significant after FDR correction.
- High-beta contrast after removal of the strongest suppressor: −0.943, exact p = 0.06838.
- Fractional high-beta contrast (20–30 Hz / 0.5–30 Hz): −0.02268, exact p = 0.14188.
- PCA explained variance: PC1 74.129%, PC2 14.875%, PC3 8.875%, PC4 2.121%.
- Four-feature multivariate group comparison: exact PERMANOVA p = 0.42595.

## Reproducibility and access

The repository includes the 300-Hz analysis script, aggregate-only results, methodological documentation, and synthetic numerical tests. Replication from source recordings requires authorized access to the original study data. Raw EEG, participant-level derived records, BrainAvatar or cognitive-test workbooks, participant crosswalks, and individual PCA coordinates are excluded from this release.

Tag: `v2.0.0`. The GitHub–Zenodo integration should archive this release as a new version and assign its own version-specific DOI.
