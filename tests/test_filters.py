import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from cellqc import compute_drops, filter_cells, flag_cells, sensitivity  # noqa: E402
from make_synthetic_data import make  # noqa: E402


def tiny():
    return pd.DataFrame({
        "condition": ["A", "A", "A", "B"],
        "r_0": [100.0, 100.0, 0.0, 100.0],
        "r_3": [90.0, 50.0, 10.0, 100.0],
        "r_5": [85.0, 45.0, 9.5, 20.0],
    })


def test_percent_drop_and_zero_baseline():
    d = compute_drops(tiny())
    assert d.loc[0, "pct_drop_r0_r3"] == pytest.approx(0.10)
    assert d.loc[3, "pct_drop_r3_r5"] == pytest.approx(0.80)
    assert np.isnan(d.loc[2, "pct_drop_r0_r3"])  # r_0 == 0 -> NaN, not inf


def test_fixed_threshold_flags_either_window():
    annotated, summary = flag_cells(tiny(), method="fixed", threshold=0.2)
    assert annotated["flag_any"].tolist() == [False, True, False, True]
    assert len(filter_cells(annotated)) == 2
    a = summary.set_index("condition")
    assert a.loc["A", "n_flagged_any"] == 1
    assert a.loc["B", "fraction_removed"] == pytest.approx(1.0)


def test_ksd_threshold_is_mean_plus_k_sd_per_group():
    df = make(200, seed=1)
    _, summary = flag_cells(df, method="ksd", k=2)
    for _, row in summary.iterrows():
        assert row["threshold_r0_r3"] == pytest.approx(row["mean_r0_r3"] + 2 * row["sd_r0_r3"])


def test_recovers_synthetic_ground_truth():
    df = make(500, seed=2)
    annotated, _ = flag_cells(df, method="fixed", threshold=0.2)
    truth, pred = annotated["true_bad"], annotated["flag_any"]
    recall = (truth & pred).sum() / truth.sum()
    precision = (truth & pred).sum() / pred.sum()
    assert recall > 0.95 and precision > 0.9


def test_sensitivity_removal_decreases_with_k():
    df = make(300, seed=3)
    s = sensitivity(df, values=[1, 2, 3]).groupby("parameter")["fraction_removed"].mean()
    assert s.is_monotonic_decreasing
