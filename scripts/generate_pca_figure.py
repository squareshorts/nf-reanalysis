#!/usr/bin/env python3
"""Generate the four-feature PCA figure for YNIRP-D-26-00004.

The script refuses to run unless subject_level_eeg.csv matches the SHA256
reported by the completed local revision analysis. This prevents a synthetic
or different participant file from being used in the manuscript figure.
"""
from __future__ import annotations
import argparse
import hashlib
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

EXPECTED_SHA256 = "7c2d695b51a152fe1036e7a8d60bb8b541301a2700c1a04a97adab49cc7664c0"
EXPECTED_VARIANCE = np.array([69.9831, 20.4963, 9.1691, 0.3514])
FEATURES = ["theta_power", "smr_power", "highbeta_power", "pac_z"]
REFERENCE_LOADINGS = np.array([
    [0.5549, -0.1521, -0.5576, -0.5984],
    [0.5833, -0.1752, -0.1964,  0.7684],
    [0.5103, -0.1986,  0.8054, -0.2268],
    [0.3024,  0.9522,  0.0428, -0.0014],
])


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Authorized local subject-level derived EEG CSV")
    ap.add_argument("--output", required=True, help="Local PNG output path; review confidentiality before sharing")
    args = ap.parse_args()

    data_path = Path(args.input).expanduser().resolve()
    if not data_path.is_file():
        raise FileNotFoundError(f"Missing input: {data_path}")

    actual_hash = sha256(data_path)
    if actual_hash != EXPECTED_SHA256:
        raise RuntimeError(
            "Input SHA256 does not match the dataset used for the completed "
            "revision report. Figure generation stopped.\n"
            f"Expected: {EXPECTED_SHA256}\nActual:   {actual_hash}"
        )

    df = pd.read_csv(data_path)
    required = {"subject", "group"}
    for feature in FEATURES:
        required.update({f"{feature}_pre", f"{feature}_post"})
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    counts = df["group"].value_counts().to_dict()
    if len(df) != 16 or counts.get("trained") != 8 or counts.get("control") != 8:
        raise ValueError(f"Unexpected sample structure: rows={len(df)}, groups={counts}")

    for feature in FEATURES:
        df[f"delta_{feature}"] = df[f"{feature}_post"] - df[f"{feature}_pre"]

    trained = df.loc[df.group.eq("trained")].copy()
    control = df.loc[df.group.eq("control")].copy()
    columns = [f"delta_{feature}" for feature in FEATURES]
    X_trained = trained[columns].to_numpy(dtype=float)
    X_control = control[columns].to_numpy(dtype=float)

    scaler = StandardScaler().fit(X_trained)
    Z_trained = scaler.transform(X_trained)
    pca = PCA(n_components=4).fit(Z_trained)
    variance = pca.explained_variance_ratio_ * 100
    if not np.allclose(variance, EXPECTED_VARIANCE, atol=0.02, rtol=0):
        raise RuntimeError(
            "PCA variance does not reproduce the audited four-feature solution. "
            f"Observed={np.round(variance, 4).tolist()}"
        )

    # Align arbitrary PCA component signs to the audited loading table.
    for j in range(4):
        if np.dot(pca.components_[j], REFERENCE_LOADINGS[:, j]) < 0:
            pca.components_[j] *= -1
    scores = pca.transform(Z_trained)

    delta_hb_control = control["delta_highbeta_power"].to_numpy(dtype=float)
    delta_hb_trained = trained["delta_highbeta_power"].to_numpy(dtype=float)
    ref_mean = float(delta_hb_control.mean())
    ref_sd = float(delta_hb_control.std(ddof=1))
    z_hb = (delta_hb_trained - ref_mean) / ref_sd

    categories = []
    for value in z_hb:
        if value <= -2:
            categories.append("Strong suppressor")
        elif value <= -1:
            categories.append("Moderate suppressor")
        elif value < 1:
            categories.append("Non-responder")
        else:
            categories.append("High-beta increaser")

    palette = {
        "Strong suppressor": "#D55E00",
        "Moderate suppressor": "#0072B2",
        "Non-responder": "#777777",
        "High-beta increaser": "#009E73",
    }
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    for category in ["Strong suppressor", "Moderate suppressor", "Non-responder", "High-beta increaser"]:
        idx = [i for i, label in enumerate(categories) if label == category]
        if idx:
            ax.scatter(
                scores[idx, 0], scores[idx, 1], s=68, label=f"{category} (n={len(idx)})",
                color=palette[category], edgecolor="black", linewidth=0.55, zorder=3
            )
    ax.axhline(0, color="0.75", lw=0.8, zorder=1)
    ax.axvline(0, color="0.75", lw=0.8, zorder=1)
    ax.set_xlabel(f"PC1 ({variance[0]:.2f}% variance)")
    ax.set_ylabel(f"PC2 ({variance[1]:.2f}% variance)")
    ax.set_title("Four-feature PCA of trained-group resting EEG changes")
    ax.legend(frameon=False, loc="best")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("PCA figure generated successfully.")
    print(f"Input SHA256: {actual_hash}")
    print(f"Variance explained (%): {np.round(variance, 4).tolist()}")
    print(f"Responder counts: {pd.Series(categories).value_counts().to_dict()}")
    print(f"Output: {output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
