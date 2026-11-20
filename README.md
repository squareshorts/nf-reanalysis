# Resting-state EEG reanalysis following SMR neurofeedback

Version **2.0.0**. Reproducibility resources for *Responder heterogeneity and multivariate organization of resting-state EEG following SMR neurofeedback* (NeuroImage: Reports, revised manuscript).

## Analysis

The analysis uses 32 pre/post C3 recording segments (16 excavator operators, trained n=8; nonfeedback control n=8), acquired at **300 Hz**. The software recomputes spectral features by Welch integration and delta-phase--SMR-amplitude coupling using the original 200-surrogate algorithm; it then runs exact two-sided permutation tests, bootstrap confidence intervals, four-feature Benjamini--Hochberg adjustment, control-referenced responder classification, PCA bootstrap/leave-one-out tests, and Euclidean PERMANOVA and PERMDISP-style tests.

The primary features are theta (4--7 Hz), SMR (12--15 Hz), high-beta (20--30 Hz), and delta-phase--SMR-amplitude coupling. Relative high-beta power uses 0.5--30 Hz as denominator.

## Reproducibility files

- `scripts/recompute_full_300hz.py`: restricted-source reproducibility implementation.
- `results/aggregate/aggregate_300hz_public.json`: results derived from the original recordings; aggregate-only.
- `tests/`: numerical smoke tests on generated data (not empirical replication).
- `docs/`: protocol and data-governance documentation.

Install with `python -m pip install -r requirements.txt`. Run numerical tests using `python -m pytest`.

To reproduce empirical results, an authorized local copy of the source EEG segments and the original participant-level feature table is required. The analysis refuses mismatched original inputs.

Example on a private computer with those data:

```powershell
python scripts/recompute_full_300hz.py --root "C:\path\to\authorized\nf_reanalysis" --outdir "C:\work\nf-reanalysis-300hz-revision"
```

The example paths are local only and can be replaced according to the user's authorized input location. Do not commit any private inputs or participant-level outputs.

## Main reported results

Absolute high-beta power: trained-minus-control change -1.452 (exact p=0.0340, four-feature FDR q=0.1361). Removing the strongest suppressor gives -0.943 (p=0.0684). Relative high-beta change: -0.02268 (p=0.1419). Four-feature PERMANOVA p=0.4260. PCA variance: 74.129%, 14.875%, 8.875%, 2.121%.

The distinctions between nominal, multiplicity-adjusted and exploratory tests are essential to interpretation. The nonfeedback control does not establish the specificity of contingent reinforcement.

## Availability

Original EDF/MAT recordings, individual-level EEG feature tables, BrainAvatar scoring files, cognitive assessments, source-record crosswalks and individual PCA coordinates are excluded from this repository due to institutional confidentiality restrictions. Data access may be considered through the corresponding author under institutional approval.

## Citation

Cite the version-specific Zenodo DOI for the **v2.0.0** snapshot once released. The earlier `v1.0.0` DOI `10.5281/zenodo.24378887` refers to the superseded 256-Hz implementation and should not be used to support the 300-Hz results.

## Authorship

Suzana Cescon de Souza; Caio de Oliveira Bastos; Arthur Gonsales da Silva; Sergio Ivan Viademonte; Ana Carolina Quintão Siravenha; Bruno Duarte Gomes; Antonio Pereira Jr.

Correspondence: apereira@ufpa.br

## License

A separate software reuse license has not been assigned. Public access to the archive does not imply unrestricted software reuse.
