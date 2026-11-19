#!/usr/bin/env python3
"""Aggregate statistical reanalysis for the occupational SMR neurofeedback EEG study.

The script consumes an authorized local subject-level derived EEG table. It does
not include or process raw EDF/MAT files, BrainAvatar workbooks, or participant
crosswalks. It writes aggregate tables and a local-only unlabeled PCA plot.

Do not commit the input table or local run directory.
"""
from __future__ import annotations

import argparse
import itertools
import hashlib
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import linear_sum_assignment
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

SEED = 20261117
FEATURES = ("theta_power", "smr_power", "highbeta_power", "pac_z")
DISPLAY = {"theta_power": "Theta (4-7 Hz)", "smr_power": "SMR (12-15 Hz)",
           "highbeta_power": "High-beta (20-30 Hz)", "pac_z": "Delta-SMR PAC (z)"}
N_EXPECTED = 8
MIN_SUBGROUP_FOR_AGGREGATES = 3


def say(lines, message=""):
    lines.append(str(message))


def validate_data(path):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Missing input: {path}")
    df = pd.read_csv(path)
    required = {"subject", "group"} | {f"{f}_{t}" for f in FEATURES for t in ("pre", "post")}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Required columns missing from the supplied derived EEG table: {missing}")
    df["group"] = df["group"].astype(str).str.strip().str.lower()
    if set(df["group"]) != {"control", "trained"}:
        raise ValueError("Unexpected group labels; expected trained and control")
    counts = df["group"].value_counts().to_dict()
    if counts.get("control") != N_EXPECTED or counts.get("trained") != N_EXPECTED:
        raise ValueError(f"Expected 8 participants per group; found {counts}")
    if df["subject"].isna().any() or df["subject"].duplicated().any():
        raise ValueError("Missing or duplicated subject identifiers")
    for f in FEATURES:
        for t in ("pre", "post"):
            col = f"{f}_{t}"
            df[col] = pd.to_numeric(df[col], errors="coerce")
            if not np.isfinite(df[col]).all():
                raise ValueError(f"Missing / nonnumeric / nonfinite EEG values in {col}")
            if f != "pac_z" and (df[col] < 0).any():
                raise ValueError(f"Negative spectral power in {col}")
        df[f"change_{f}"] = df[f"{f}_post"] - df[f"{f}_pre"]
    return df


def exact_assignment_indices(n, n_tr):
    return np.array(list(itertools.combinations(range(n), n_tr)), dtype=np.int16)


def exact_did_p(x, observed_tr, assignments):
    n = len(x)
    n_tr = observed_tr.sum()
    k = len(assignments)
    # Group means, centered by pooled mean; membership varies across all allocations.
    mask = np.zeros((k, n), dtype=np.int8)
    mask[np.arange(k)[:, None], assignments] = 1
    sum_tr = mask @ x
    means = sum_tr / n_tr - (np.sum(x) - sum_tr) / (n - n_tr)
    obs = x[observed_tr].mean() - x[~observed_tr].mean()
    p = np.mean(np.abs(means) >= abs(obs) - 1e-12)
    return float(obs), float(p)


def bh_qvalues(p):
    p = np.asarray(p, float)
    ix = np.argsort(p)
    qsort = np.minimum.accumulate((p[ix] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(p)
    out[ix] = np.clip(qsort, 0, 1)
    return out


def group_and_interactions(df, outdir, lines, rng, assignments, n_boot):
    rows, prepost = [], []
    for f in FEATURES:
        for g in ("control", "trained"):
            sub = df[df.group == g]
            pre, post = sub[f"{f}_pre"].to_numpy(), sub[f"{f}_post"].to_numpy()
            diffs = post - pre
            for t, vals in (("pre", pre), ("post", post)):
                prepost.append({"feature": f, "group": g, "time": t, "n": len(vals),
                                "mean": float(vals.mean()), "sd": float(vals.std(ddof=1))})
        tr = df.loc[df.group == "trained", f"change_{f}"].to_numpy()
        co = df.loc[df.group == "control", f"change_{f}"].to_numpy()
        combined = df[f"change_{f}"].to_numpy()
        observed_tr = (df.group == "trained").to_numpy()
        contrast, pexact = exact_did_p(combined, observed_tr, assignments)
        wt = stats.ttest_ind(tr, co, equal_var=False)
        mw = stats.mannwhitneyu(tr, co, alternative="two-sided")
        pooled_sd = math.sqrt(((len(tr) - 1) * tr.var(ddof=1) + (len(co) - 1) * co.var(ddof=1)) /
                              (len(tr) + len(co) - 2))
        d = contrast / pooled_sd if pooled_sd > 0 else np.nan
        dfree = len(tr) + len(co) - 2
        hedges_g = d * (1 - 3 / (4 * dfree - 1))
        group_indices_tr = rng.integers(0, len(tr), size=(n_boot, len(tr)))
        group_indices_co = rng.integers(0, len(co), size=(n_boot, len(co)))
        ci_boot = tr[group_indices_tr].mean(axis=1) - co[group_indices_co].mean(axis=1)
        ci_lo, ci_hi = np.percentile(ci_boot, (2.5, 97.5))
        # Two-timepoint group x time coefficient equals the between-group difference in changes.
        rows.append({"feature": f, "trained_change_mean": tr.mean(),
                     "control_change_mean": co.mean(),
                     "group_x_time_difference_of_changes": contrast,
                     "bootstrap_95ci_low": ci_lo, "bootstrap_95ci_high": ci_hi,
                     "exact_assignment_two_sided_p": pexact,
                     "welch_change_t": float(wt.statistic), "welch_change_p": float(wt.pvalue),
                     "mann_whitney_p": float(mw.pvalue), "change_hedges_g": hedges_g,
                     "change_rank_biserial": 2 * float(mw.statistic) / (len(tr) * len(co)) - 1})
    tab = pd.DataFrame(rows)
    tab["exact_p_BH_FDR_4_features"] = bh_qvalues(tab.exact_assignment_two_sided_p)
    tab.to_csv(outdir / "group_time_interaction.csv", index=False)
    pd.DataFrame(prepost).to_csv(outdir / "prepost_aggregate.csv", index=False)
    say(lines, "GROUP x TIME: difference of post-pre changes (trained minus control)")
    say(lines, "Exact two-sided allocation p-values enumerate 12,870 group assignments; "
               "95% CIs are stratified subject-level bootstrap intervals.")
    say(lines, tab.round(5).to_string(index=False))
    say(lines, "Interpretation: this is the 2-timepoint interaction contrast; no independence "
               "assumption is made across the repeated pre/post observations.")
    say(lines, "FDR adjustment treats all four EEG features as a family; identify any "
               "prospectively designated primary endpoint separately in the manuscript.")
    say(lines)
    return tab


def category(z, strong=-2.0, moderate=-1.0):
    if z <= strong:
        return "strong_suppressor"
    if z <= moderate:
        return "moderate_suppressor"
    if z < 1.0:
        return "non_responder"
    return "highbeta_increaser"


def responder_analysis(df, outdir, lines):
    dctrl = df.loc[df.group == "control", "change_highbeta_power"].to_numpy()
    mu, sd = float(np.mean(dctrl)), float(np.std(dctrl, ddof=1))
    if sd <= 0 or not np.isfinite(sd):
        raise ValueError("Control high-beta changes have no valid reference SD")
    tr = df.loc[df.group == "trained"].copy()
    tr["hb_z"] = (tr.change_highbeta_power - mu) / sd
    base = [category(z) for z in tr.hb_z]
    categories = ("strong_suppressor", "moderate_suppressor", "non_responder", "highbeta_increaser")
    rows = []
    for strong in (-2.5, -2.0, -1.5):
        for moderate in (-1.5, -1.0, -0.5):
            if strong >= moderate:
                continue
            current = [category(z, strong, moderate) for z in tr.hb_z]
            row = {"z_strong_cutoff": strong, "z_moderate_cutoff": moderate,
                   "n_changed_vs_original": sum(a != b for a, b in zip(current, base))}
            row.update({"n_" + cat: current.count(cat) for cat in categories})
            rows.append(row)
    pd.DataFrame(rows).to_csv(outdir / "responder_threshold_sensitivity.csv", index=False)
    tr["category"] = base
    baseline_rows = []
    for cat in categories:
        sub = tr.loc[tr.category == cat]
        if len(sub) == 0:
            continue
        vals = sub.highbeta_power_pre.to_numpy()
        row = {"responder_category": cat, "n": len(vals),
               "baseline_hb_mean": float(vals.mean()) if len(vals) >= MIN_SUBGROUP_FOR_AGGREGATES else np.nan,
               "baseline_hb_sd": float(vals.std(ddof=1)) if len(vals) >= MIN_SUBGROUP_FOR_AGGREGATES else np.nan,
               "suppressed_for_small_group": len(vals) < MIN_SUBGROUP_FOR_AGGREGATES}
        baseline_rows.append(row)
    pd.DataFrame(baseline_rows).to_csv(outdir / "responder_baseline_aggregate.csv", index=False)
    say(lines, f"RESPONDERS: control high-beta change reference mean={mu:.5f}, sample SD={sd:.5f}")
    say(lines, "Primary counts: " + ", ".join(f"{cat}={base.count(cat)}" for cat in categories))
    say(lines, "Threshold shifts of +/-0.5 SD (including combinations) saved to CSV.")
    say(lines, "Baseline high-beta subgroup means are suppressed for n < 3 to reduce disclosure risk.")
    say(lines, "Age was not available in the supplied derived EEG table; no responder-subgroup age comparisons were calculated.")
    say(lines)
    return tr


def align_components(reference_components, new_components):
    similarity = reference_components @ new_components.T
    ref_idx, new_idx = linear_sum_assignment(-np.abs(similarity))
    components = np.empty_like(new_components)
    matched = np.zeros(reference_components.shape[0], dtype=int)
    cos = np.zeros(reference_components.shape[0], dtype=float)
    for i, j in zip(ref_idx, new_idx):
        sign = 1.0 if similarity[i, j] >= 0 else -1.0
        components[i] = sign * new_components[j]
        matched[i] = j
        cos[i] = abs(similarity[i, j])
    return components, matched, cos


def fit_pca(x):
    s = StandardScaler()
    z = s.fit_transform(x)
    if np.any(s.scale_ <= 1e-12):
        raise ValueError("Zero-variance predictor in PCA")
    pca = PCA(n_components=len(FEATURES)).fit(z)
    return pca, z


def pca_analysis(df, tr, outdir, lines, rng, n_boot):
    names = ["change_" + f for f in FEATURES]
    X = tr[names].to_numpy()
    pca, z = fit_pca(X)
    loadings = pd.DataFrame(pca.components_.T, index=FEATURES,
                            columns=[f"PC{i}" for i in range(1, 5)])
    loadings.index.name = "feature"
    loadings.to_csv(outdir / "pca_four_feature_loadings.csv")
    ev = pd.DataFrame({"component": [f"PC{i}" for i in range(1, 5)],
                       "variance_percent": 100 * pca.explained_variance_ratio_,
                       "cumulative_percent": 100 * np.cumsum(pca.explained_variance_ratio_)})
    ev.to_csv(outdir / "pca_four_feature_variance.csv", index=False)
    say(lines, "PCA: trained-group within-sample standardized four-feature changes")
    say(lines, ev.round(4).to_string(index=False))
    say(lines, loadings.round(4).to_string())

    # Bootstrap refits the scaling and all four components on every training-group resample.
    boot_components = []
    boot_cosines = []
    boot_ev = []
    skipped = 0
    for _ in range(n_boot):
        idx = rng.integers(0, len(X), size=len(X))
        subset = X[idx]
        if np.any(np.std(subset, axis=0, ddof=0) <= 1e-12):
            skipped += 1
            continue
        p, _ = fit_pca(subset)
        aligned, assignment, cos = align_components(pca.components_, p.components_)
        boot_components.append(aligned[:2].T)
        boot_cosines.append(cos[:2])
        boot_ev.append(p.explained_variance_ratio_[assignment[:2]])
    if not boot_components:
        raise RuntimeError("No valid PCA bootstrap samples")
    boots = np.stack(boot_components)
    cosine = np.array(boot_cosines)
    boot_var = np.asarray(boot_ev)
    boot_rows = []
    for pc in range(2):
        for j, f in enumerate(FEATURES):
            arr = boots[:, j, pc]
            boot_rows.append({"component": f"PC{pc+1}", "feature": f,
                              "point_loading": float(pca.components_[pc, j]),
                              "bootstrap_2.5pct": float(np.percentile(arr, 2.5)),
                              "bootstrap_median": float(np.median(arr)),
                              "bootstrap_97.5pct": float(np.percentile(arr, 97.5))})
    pd.DataFrame(boot_rows).to_csv(outdir / "pca_bootstrap_loadings.csv", index=False)

    loo_rows = []
    for i in range(len(X)):
        reduced = np.delete(X, i, axis=0)
        p, _ = fit_pca(reduced)
        _, assignment, cos = align_components(pca.components_, p.components_)
        loo_rows.append({"pc1_cosine": cos[0], "pc2_cosine": cos[1],
                         "pc1_variance_percent": 100*p.explained_variance_ratio_[assignment[0]],
                         "pc2_variance_percent": 100*p.explained_variance_ratio_[assignment[1]]})
    loo = pd.DataFrame(loo_rows)
    # Never write an identifiable leave-one-out row. Aggregate ranges and medians only.
    summaries = []
    for key in loo.columns:
        summaries.append({"metric": key, "min": float(loo[key].min()),
                          "median": float(loo[key].median()), "max": float(loo[key].max())})
    pd.DataFrame(summaries).to_csv(outdir / "pca_leave_one_out_aggregate.csv", index=False)
    say(lines, f"Bootstrap PCA refits={len(boots)} (discarded degenerate resamples={skipped})")
    say(lines, "Bootstrap aligned-loading interval table saved to pca_bootstrap_loadings.csv.")
    for pc in (0, 1):
        say(lines, f"PC{pc+1} cosine similarity: bootstrap median={np.median(cosine[:,pc]):.3f}, "
                   f"95% percentile interval={np.percentile(cosine[:,pc],[2.5,97.5]).round(3).tolist()}, "
                   f"LOOCV range={loo[f'pc{pc+1}_cosine'].min():.3f}..{loo[f'pc{pc+1}_cosine'].max():.3f}")
        say(lines, f"PC{pc+1} aligned bootstrap variance percentiles: "
                   f"{np.percentile(100*boot_var[:,pc],[2.5,50,97.5]).round(3).tolist()}")

    # Sensitivity only: using reference scaling from the control group changes PCA geometry.
    Xcontrol = df.loc[df.group == "control", names].to_numpy()
    ctrl_mu, ctrl_sd = Xcontrol.mean(axis=0), Xcontrol.std(axis=0, ddof=1)
    if np.any(ctrl_sd <= 1e-12):
        raise ValueError("Control-reference PCA scaling includes zero SD")
    alt = PCA(n_components=4).fit((X - ctrl_mu) / ctrl_sd)
    pd.DataFrame({"component": [f"PC{i}" for i in range(1, 5)],
                  "control_ref_variance_percent": 100*alt.explained_variance_ratio_}).to_csv(
        outdir / "pca_control_reference_sensitivity.csv", index=False)
    say(lines, "Control-referenced FOUR-feature PCA sensitivity variance percentages: " +
               str(np.round(100 * alt.explained_variance_ratio_, 3).tolist()))
    say(lines)

    # Local-only visualization: unlabeled points, no subject IDs, unlike original figure.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    score = pca.transform(z)
    colors = {"strong_suppressor": "#b2182b", "moderate_suppressor": "#ef8a62",
              "non_responder": "#888888", "highbeta_increaser": "#2166ac"}
    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    for k, cat in enumerate(tr["category"].to_numpy()):
        ax.scatter(score[k, 0], score[k, 1], color=colors[cat], s=65, edgecolors="black", linewidths=0.5)
    for lab in dict.fromkeys(tr["category"].tolist()):
        ax.scatter([], [], s=65, color=colors[lab], label=lab.replace("_", " "))
    ax.axhline(0, color="gray", lw=0.7)
    ax.axvline(0, color="gray", lw=0.7)
    ax.set(xlabel=f"PC1 ({100*pca.explained_variance_ratio_[0]:.1f}% variance)",
           ylabel=f"PC2 ({100*pca.explained_variance_ratio_[1]:.1f}% variance)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "figures_local_only" / "fig_pca_four_feature_corrected.png", dpi=300)
    plt.close(fig)
    return pd.DataFrame({"subject": tr.subject.to_numpy(), "PC1": score[:, 0],
                         "hb_z": tr.hb_z.to_numpy(), "delta_hb": tr.change_highbeta_power.to_numpy()})


def permanova(df, outdir, lines, assignments):
    # Pooled (group-label-independent) scaling is necessary for exact label permutation.
    x = df[["change_" + f for f in FEATURES]].to_numpy()
    x = StandardScaler().fit_transform(x)
    n, nvar = x.shape
    is_tr = (df.group == "trained").to_numpy()
    n_tr = int(is_tr.sum())
    n_ctrl = n - n_tr
    ss_total = float(np.sum((x - x.mean(axis=0))**2))
    mask = np.zeros((len(assignments), n), dtype=np.int8)
    mask[np.arange(len(assignments))[:, None], assignments] = 1
    means_tr = (mask @ x) / n_tr
    means_co = (np.sum(x, axis=0) - mask @ x) / n_ctrl
    ss_between = n_tr * n_ctrl / n * np.sum((means_tr - means_co)**2, axis=1)
    ss_within = ss_total - ss_between
    f_perm = (ss_between / 1) / (ss_within / (n - 2))
    center_diff = x[is_tr].mean(axis=0) - x[~is_tr].mean(axis=0)
    obs_between = n_tr * n_ctrl / n * np.dot(center_diff, center_diff)
    obs_within = ss_total - obs_between
    obs_F = obs_between / (obs_within / (n - 2))
    pexact = float(np.mean(f_perm >= obs_F - 1e-12))
    r2 = obs_between / ss_total
    # Diagnostic: group dispersion difference with group centroids refit for each allocation.
    # This is not a formal PERMDISP F-test; descriptive exact-label assessment only.
    def dispersion_difference(tr_mask):
        mt = x[tr_mask].mean(axis=0)
        mc = x[~tr_mask].mean(axis=0)
        dt = np.linalg.norm(x[tr_mask] - mt, axis=1).mean()
        dc = np.linalg.norm(x[~tr_mask] - mc, axis=1).mean()
        return float(dt-dc), float(dt), float(dc)
    d_obs, dt_obs, dc_obs = dispersion_difference(is_tr)
    d_perm = np.zeros(len(assignments))
    for k, idx in enumerate(assignments):
        bm = np.zeros(n, dtype=bool)
        bm[idx] = True
        d_perm[k] = dispersion_difference(bm)[0]
    pdisp = float(np.mean(np.abs(d_perm) >= abs(d_obs) - 1e-12))
    out = {"distance": "pooled-standardized Euclidean (4 EEG change features)",
           "n": n, "n_trained": n_tr, "n_control": n_ctrl,
           "pseudo_F": obs_F, "exact_allocation_p": pexact, "R2": r2,
           "n_assignments": len(assignments),
           "mean_distance_centroid_trained": dt_obs,
           "mean_distance_centroid_control": dc_obs,
           "dispersion_mean_difference": d_obs,
           "exploratory_dispersion_exact_p": pdisp}
    pd.DataFrame([out]).to_csv(outdir / "permanova_corrected.csv", index=False)
    say(lines, "FOUR-FEATURE PERMANOVA: pooled-standardized Euclidean distances; exact allocations")
    say(lines, f"pseudo-F={obs_F:.6f}, R2={r2:.6f}, exact p={pexact:.6f}, assignments={len(assignments)}")
    say(lines, f"Exploratory centroid-dispersion difference={d_obs:.6f}, exact label p={pdisp:.6f}")
    say(lines, "PERMANOVA is sensitive to distributional dispersion; do not interpret the p-value "
               "as exclusively a centroid test when dispersions differ.")
    say(lines)
    return out





def exact_permdisp_euclidean(df, outdir, lines, assignments):
    """Exact two-group dispersion test using distances to group centroids.

    Distances are computed in pooled-standardized four-feature change space.
    For every balanced allocation, centroids and distances are recomputed and
    the one-way ANOVA F statistic on distances is recalculated.
    """
    cols = ["change_" + f for f in FEATURES]
    x = StandardScaler().fit_transform(df[cols].to_numpy(dtype=float))
    labels = (df.group == "trained").to_numpy()
    n = len(df)

    def statistic(mask):
        groups = [x[mask], x[~mask]]
        distances = [np.linalg.norm(g - g.mean(axis=0), axis=1) for g in groups]
        all_d = np.concatenate(distances)
        grand = all_d.mean()
        ss_between = sum(len(d) * (d.mean() - grand) ** 2 for d in distances)
        ss_within = sum(np.sum((d - d.mean()) ** 2) for d in distances)
        f_value = (ss_between / 1) / (ss_within / (n - 2)) if ss_within > 0 else np.inf
        return float(f_value), float(distances[0].mean()), float(distances[1].mean())

    observed_f, mean_trained, mean_control = statistic(labels)
    permuted_f = np.empty(len(assignments), dtype=float)
    for k, idx in enumerate(assignments):
        mask = np.zeros(n, dtype=bool)
        mask[idx] = True
        permuted_f[k] = statistic(mask)[0]
    p_exact = float(np.mean(permuted_f >= observed_f - 1e-12))
    result = {
        "distance": "pooled-standardized Euclidean (four EEG change features)",
        "n": n,
        "n_trained": int(labels.sum()),
        "n_control": int((~labels).sum()),
        "permdisp_F": observed_f,
        "exact_upper_tail_p": p_exact,
        "mean_distance_trained": mean_trained,
        "mean_distance_control": mean_control,
        "n_assignments": len(assignments),
    }
    pd.DataFrame([result]).to_csv(outdir / "permdisp_exact.csv", index=False)
    say(lines, "FORMAL EUCLIDEAN PERMDISP-STYLE TEST (distances to group centroids)")
    say(lines, f"Pooled-standardized four-feature changes; F={observed_f:.6f}; exact upper-tail allocation p={p_exact:.6f}; assignments={len(assignments)}.")
    say(lines, f"Mean distance to centroid: trained={mean_trained:.6f}; control={mean_control:.6f}; difference={mean_trained-mean_control:.6f}.")
    say(lines, "This exact upper-tail allocation test is reported separately from the PERMANOVA group-centroid test.")
    say(lines)
    return result


def bayes_optional(df, outdir, lines):
    try:
        import pingouin as pg
    except ImportError:
        say(lines, "BAYES FACTORS: not recomputed (optional package pingouin not installed).")
        say(lines, "Submitted BF01 values are reciprocals of submitted BF10 values: "
                   "SMR: 2.00 (BF10 0.50); PAC: 1.33 (BF10 0.75). "
                   "These are NOT new statistical results.")
        return
    tr = df[df.group == "trained"]
    rows = []
    for f in ("smr_power", "pac_z"):
        res = pg.ttest(tr[f"{f}_post"], tr[f"{f}_pre"], paired=True)
        bf10 = float(res.BF10.iloc[0])
        rows.append({"feature": f, "bf10": bf10, "bf01": 1/bf10})
    pd.DataFrame(rows).to_csv(outdir / "bayes_factors.csv", index=False)
    say(lines, "BAYES FACTORS: recomputed using locally installed Pingouin defaults")
    say(lines, pd.DataFrame(rows).round(5).to_string(index=False))
    say(lines)



def main():
    ap = argparse.ArgumentParser(description="Reproduce aggregate analyses for the SMR neurofeedback EEG reanalysis")
    ap.add_argument("--input", type=Path, required=True,
                    help="Authorized local subject-level derived EEG CSV; never commit this file")
    ap.add_argument("--outdir", type=Path, required=True,
                    help="Local output directory; should remain outside the public repository")
    ap.add_argument("--boots", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--expected-sha256", type=str, default=None,
                    help="Optional expected SHA256 for strict input provenance validation")
    ap.add_argument("--skip-bayes", action="store_true",
                    help="Skip paired-sample Bayes factors if Pingouin is unavailable")
    args = ap.parse_args()
    if args.boots < 200:
        ap.error("--boots must be at least 200")

    input_path = args.input.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    input_hash = hashlib.sha256(input_path.read_bytes()).hexdigest()
    if args.expected_sha256 and input_hash.lower() != args.expected_sha256.strip().lower():
        raise ValueError("Input SHA256 does not match the expected source dataset. "
                         f"Expected={args.expected_sha256}; observed={input_hash}")

    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "figures_local_only").mkdir(exist_ok=True)
    owned_files = [
        "revision_aggregate_report.txt", "group_time_interaction.csv",
        "prepost_aggregate.csv", "responder_threshold_sensitivity.csv",
        "responder_baseline_aggregate.csv", "pca_four_feature_loadings.csv",
        "pca_four_feature_variance.csv", "pca_bootstrap_loadings.csv",
        "pca_leave_one_out_aggregate.csv", "pca_control_reference_sensitivity.csv",
        "permanova_corrected.csv", "permdisp_exact.csv", "bayes_factors.csv",
    ]
    for name in owned_files:
        (outdir / name).unlink(missing_ok=True)
    (outdir / "figures_local_only" / "fig_pca_four_feature_corrected.png").unlink(missing_ok=True)

    rng = np.random.default_rng(args.seed)
    lines = []
    say(lines, "SMR NEUROFEEDBACK EEG REANALYSIS - AGGREGATE REVISION REPORT")
    say(lines, "Outputs contain aggregate statistics only; no participant identifiers or participant-level tables are written.")
    say(lines, f"Input SHA256: {input_hash}")
    say(lines, f"Reproducibility: seed={args.seed}; bootstrap iterations={args.boots}")
    say(lines)
    df = validate_data(input_path)
    n_tr = int((df.group == "trained").sum())
    assignments = exact_assignment_indices(len(df), n_tr)
    say(lines, f"Validated participants: trained={n_tr}, control={len(df)-n_tr}; exact balanced allocations={len(assignments)}")
    say(lines)

    group_and_interactions(df, outdir, lines, rng, assignments, args.boots)
    trained = responder_analysis(df, outdir, lines)
    pca_analysis(df, trained, outdir, lines, rng, args.boots)
    permanova(df, outdir, lines, assignments)
    exact_permdisp_euclidean(df, outdir, lines, assignments)
    if args.skip_bayes:
        say(lines, "BAYES FACTORS: skipped by user request.")
        say(lines)
    else:
        bayes_optional(df, outdir, lines)
    say(lines, "ONLINE/OFFLINE ASSOCIATIONS: not analyzed. Available BrainAvatar score records were not source-verified against EEG participant records; no correlations are reported.")
    say(lines, "The analysis consumes a derived participant-level EEG table. Raw EEG preprocessing and feature extraction are outside this release pipeline.")
    report = outdir / "revision_aggregate_report.txt"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Aggregate report: {report}")
    print(f"Aggregate tables: {outdir}")
    print("The PCA score plot is in figures_local_only/ and must not be published without confidentiality review.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
