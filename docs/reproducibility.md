# Reproducibility workflow

1. Install the Python dependencies from `requirements.txt`.
2. Obtain the derived EEG table through an authorized data-governance route.
3. Verify its SHA256 against the hash recorded in the manuscript revision report, if reproducing the exact archived analysis.
4. Run `scripts/run_revision_analyses.py` with an explicit input path and an output directory outside the repository.
5. Inspect `revision_aggregate_report.txt` and the aggregate CSV outputs.
6. Run `python -m pytest` to verify the pipeline on synthetic data.

The synthetic fixture validates software execution and output structure only; it does not reproduce the empirical findings. The public repository contains aggregate outputs from the empirical analysis but not the participant-level input table.


The derived EEG input used for the archived aggregate results had SHA256 `7c2d695b51a152fe1036e7a8d60bb8b541301a2700c1a04a97adab49cc7664c0`. This hash is provided solely as a file-integrity check; the underlying table is not distributed.
