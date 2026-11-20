# Reproducibility

The empirical 300-Hz analysis requires the private source EEG MAT segments, the original subject-level EEG table and the precomputed source feature table to verify input correspondence. `scripts/recompute_full_300hz.py` performs consistency checks and outputs aggregate statistical summaries.

The code's standard test suite uses synthetic signals and does not reproduce empirical results. Public aggregate results are in `results/aggregate/aggregate_300hz_public.json`; confidential observations are not distributed.
