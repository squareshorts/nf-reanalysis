from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def test_synthetic_pipeline_writes_aggregate_outputs(tmp_path):
    fixture = ROOT / "tests" / "fixtures" / "synthetic_subject_level_eeg.csv"
    out = tmp_path / "run"
    cmd = [
        sys.executable,
        str(ROOT / "scripts" / "run_revision_analyses.py"),
        "--input", str(fixture),
        "--outdir", str(out),
        "--boots", "200",
        "--skip-bayes",
    ]
    completed = subprocess.run(cmd, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
    expected = [
        "revision_aggregate_report.txt",
        "group_time_interaction.csv",
        "responder_threshold_sensitivity.csv",
        "responder_baseline_aggregate.csv",
        "pca_four_feature_loadings.csv",
        "pca_four_feature_variance.csv",
        "pca_bootstrap_loadings.csv",
        "pca_leave_one_out_aggregate.csv",
        "pca_control_reference_sensitivity.csv",
        "permanova_corrected.csv",
        "permdisp_exact.csv",
    ]
    for name in expected:
        assert (out / name).is_file(), name
    assert (out / "figures_local_only" / "fig_pca_four_feature_corrected.png").is_file()
    report = (out / "revision_aggregate_report.txt").read_text(encoding="utf-8")
    assert "no participant identifiers" in report.lower()
    assert "ONLINE/OFFLINE ASSOCIATIONS: not analyzed" in report
    # No output table may contain a participant-level identifier column.
    for csv_path in out.glob("*.csv"):
        header = csv_path.read_text(encoding="utf-8").splitlines()[0].lower()
        assert "subject" not in header
