# SMR neurofeedback EEG reanalysis

**Software release:** v1.0.0  
**Associated manuscript:** *Responder heterogeneity and multivariate organization of resting-state EEG following SMR neurofeedback* (under revision at *NeuroImage: Reports*).

This repository contains the analysis code and aggregate results for a reanalysis of resting-state EEG collected before and after an occupational sensorimotor-rhythm (SMR) neurofeedback protocol. It is a reproducibility archive for the derived-feature analyses reported in the manuscript.

## Scope

The analysis consumes a locally held participant-level derived EEG table with pre/post values for theta power, SMR power, high-beta power, and delta-phase--SMR-amplitude coupling. It does **not** recompute raw EEG preprocessing, spectral feature extraction, or PAC estimation. Those operations belong to the original acquisition/processing workflow and require access to restricted source recordings.

The public archive contains:

- the aggregate statistical analysis script;
- a guarded PCA-figure generator;
- a synthetic test fixture and automated test;
- protocol and analysis documentation;
- aggregate result tables and a sanitized aggregate report.

It does **not** contain raw EDF/MAT recordings, participant-level derived EEG rows, BrainAvatar workbooks, cognitive-test workbooks, participant crosswalks, or individual-level PCA coordinates. The PCA score plot generated from restricted data is a local output and is not included in this public archive.

## Study protocol represented in the analysis

The documented intervention comprised eight 40-minute BrainAvatar sessions over four training days, with two sessions per day separated by approximately four hours. The visual task used a spaceship and central target with an auditory reward signal. Reinforcement required increased SMR (12--15 Hz) together with reduced theta (4--7 Hz) and high-beta (20--30 Hz), sustained for at least 500 ms. The control group received no contingent feedback; no sham or yoked feedback was provided. See [`docs/protocol_summary.md`](docs/protocol_summary.md) for source-grounded detail and the documented session-count discrepancy.

## Repository layout

```text
scripts/       Reproducible analysis and local PCA figure generator
data/          Data-access instructions; no participant data
src/           Package version metadata
tests/         Synthetic-data test fixture and test
results/       Aggregate-only tables and report
docs/          Protocol, analysis, and data-governance documentation
environment/   Environment notes
```

## Environment

Python 3.10 or later is required. Create an environment and install the pinned analysis dependencies:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

On macOS/Linux, activate with `source .venv/bin/activate`.

## Reproduce the aggregate analysis

Obtain the derived EEG table through the corresponding author under the applicable institutional data-governance conditions. The required columns are `subject`, `group` (`trained` or `control`), and pre/post values for `theta_power`, `smr_power`, `highbeta_power`, and `pac_z`.

Run from the repository root, supplying the local input path and an output directory outside the tracked repository:

```powershell
python scripts/run_revision_analyses.py `
  --input "C:\path\to\authorized\subject_level_eeg.csv" `
  --outdir "C:\path\to\local_results\nf_revision_v1" `
  --boots 2000 `
  --expected-sha256 7c2d695b51a152fe1036e7a8d60bb8b541301a2700c1a04a97adab49cc7664c0
```

The SHA256 check confirms that the input is the exact derived table used for the reported revision results. If an authorized, intentionally updated input is used, omit `--expected-sha256` and document the changed input provenance. The script writes aggregate CSVs, an aggregate text report, and a local-only PCA score plot. The output directory should not be committed.

Run the synthetic-data test without access to study data:

```powershell
python -m pytest
```

The test uses artificial records and does not reproduce manuscript results.

## Principal aggregate results

- Trained-minus-control high-beta change: -1.508 (95% participant-bootstrap CI -2.720 to -0.410); exact allocation p = 0.0236; Hedges' g = -1.133; four-feature Benjamini--Hochberg q = 0.0945.
- Exact four-feature Euclidean PERMANOVA: pseudo-F = 1.1858, R2 = 0.0781, p = 0.3173.
- Exact centroid-dispersion test: F = 0.0472, p = 0.8416.
- Four-feature within-trained-group PCA variance: PC1 69.983%, PC2 20.496%, PC3 9.169%, PC4 0.351%.
- Primary high-beta categories: one strong suppressor, four moderate suppressors, and three non-responders.

All aggregate tables are in [`results/aggregate/`](results/aggregate/). Interpretation and limitations are documented in [`docs/analysis_notes.md`](docs/analysis_notes.md).

## Online performance records

Online/offline correlations are not reported. The available BrainAvatar score files did not establish a source-verified participant-level correspondence with the EEG-derived table, so no correlation estimates are included in this release.

## Data access and privacy

Raw EEG and participant-level derived data are subject to institutional and occupational confidentiality constraints and are not distributed here. BrainAvatar and cognitive-test workbooks and any participant crosswalk are also excluded. The derived EEG table may be requested from the corresponding author, subject to institutional approval and applicable data-governance requirements. See [`docs/data_governance.md`](docs/data_governance.md).

## Citation

Please cite the associated manuscript and this software release. The version-specific Zenodo DOI will be added to this README and `CITATION.cff` after the GitHub release has been archived by Zenodo.

## Authors

Suzana Cescon de Souza; Caio de Oliveira Bastos; Arthur Gonsales da Silva; Sergio Ivan Viademonte; Ana Carolina Quintão Siravenha; Bruno Duarte Gomes; Antonio Pereira Jr.

Correspondence: apereira@ufpa.br

## License

A reuse license has not yet been assigned. Until the authors add an explicit license, no permission for reuse beyond applicable statutory exceptions is granted.
