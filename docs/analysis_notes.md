# Analysis notes and interpretation

## Group-by-time contrasts

For each feature, the estimand is the trained-minus-control difference in mean post--pre change. Exact two-sided allocation p-values enumerate all 12,870 balanced assignments of eight participants to each group, consistent with the random allocation described in the source study. Stratified participant-level bootstrap intervals use 2,000 resamples. Hedges' g is reported for between-group differences in changes. Benjamini--Hochberg adjustment is applied across the four EEG features.

High-beta has the largest group-by-time contrast (difference -1.508; 95% bootstrap CI -2.720 to -0.410; exact p = 0.0236; Hedges' g = -1.133), but the four-feature adjusted q-value is 0.0945. It is a nominal contrast, not an FDR-significant result.

## Responder categories

Individual high-beta changes are standardized against the control-group change distribution. The descriptive cutoffs are z <= -2 (strong suppression), -2 < z <= -1 (moderate suppression), and -1 < z < 1 (non-responder range). These categories are not validated biological phenotypes and do not establish successful online SMR learning. Sensitivity to eight alternative cutoff combinations is reported in the aggregate tables.

## PCA

The primary PCA uses four post--pre changes, standardized within the trained group. Variance proportions are PC1 69.983%, PC2 20.496%, PC3 9.169%, and PC4 0.351%. Bootstrap loading intervals and leave-one-participant-out summaries are provided. The control-referenced scaling sensitivity changes the component ordering, so the PCA is treated as an exploratory description of covariance among the eight trained participants, not as a definitive architecture or evidence of group separation.

The PCA score plot is generated locally and is not included in this public archive because it contains eight individual data points. Review confidentiality and disclosure risk before sharing it externally.

## Multivariate group comparison

The four-feature Euclidean PERMANOVA uses pooled standardization and exact balanced allocations. The exact centroid-dispersion test separately evaluates the ANOVA F statistic on distances to group centroids. Neither test detected a group-level multivariate difference.

## Bayesian paired tests

Paired Bayes factors for trained-group SMR power and PAC are computed using Pingouin's default paired-test settings. The reported BF10/BF01 values are 0.500/2.000 for SMR and 0.745/1.342 for PAC. These values provide weak evidence and do not establish equivalence to zero.

## Online/offline associations

No online/offline correlations are included. The available BrainAvatar score records could not be linked to EEG participants through a source-verified crosswalk. No row-order or identifier-similarity matching was used.
