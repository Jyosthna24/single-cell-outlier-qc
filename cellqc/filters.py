"""Per-cell drop metrics and outlier filters for iodide-quench time courses.

During an iodide-quench experiment the red (mCherry) signal of a healthy cell
should stay roughly constant across buffer exchanges. Cells that die, detach
or are washed out of the field of view show a sudden, large red drop that has
nothing to do with CFTR activity. These functions flag such cells so they can
be excluded before CFTR function is quantified.

Two windows are checked:
    r0 -> r3  (pre-iodide / first buffer exchange)
    r3 -> r5  (after iodide addition)

A cell is flagged if its drop exceeds the threshold in EITHER window.
"""
from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import pandas as pd

WINDOWS = {
    "r0_r3": ("r_0", "r_3"),
    "r3_r5": ("r_3", "r_5"),
}


def compute_drops(df: pd.DataFrame, use_percent_drop: bool = True) -> pd.DataFrame:
    """Add absolute and fractional red-signal drops for each window.

    Fractional drop is scale-free, so it is comparable across cells with
    different expression levels. Division by zero gives NaN, not inf.
    """
    out = df.copy()
    for name, (start, end) in WINDOWS.items():
        a = out[start].astype(float)
        b = out[end].astype(float)
        out[f"delta_{name}"] = a - b
        out[f"pct_drop_{name}"] = ((a - b) / a.replace(0, np.nan)).replace(
            [np.inf, -np.inf], np.nan
        )
        out[f"metric_{name}"] = (
            out[f"pct_drop_{name}"] if use_percent_drop else out[f"delta_{name}"]
        )
    return out


def _thresholds(
    df: pd.DataFrame, group_cols: Sequence[str], method: str, threshold: float, k: float
) -> pd.DataFrame:
    """One row per group with the threshold used for each window."""
    rows = []
    for keys, g in df.groupby(list(group_cols), dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        row = dict(zip(group_cols, keys))
        for name in WINDOWS:
            v = g[f"metric_{name}"].dropna().to_numpy(float)
            n = len(v)
            mean = float(np.mean(v)) if n else float("nan")
            sd = float(np.std(v, ddof=1)) if n > 1 else 0.0
            thr = mean + k * sd if method == "ksd" else float(threshold)
            row.update({f"n_{name}": n, f"mean_{name}": mean, f"sd_{name}": sd,
                        f"threshold_{name}": thr})
        rows.append(row)
    return pd.DataFrame(rows)


def flag_cells(
    df: pd.DataFrame,
    method: str = "ksd",
    k: float = 2.0,
    threshold: float = 0.2,
    group_cols: Sequence[str] = ("condition",),
    use_percent_drop: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Flag cells with an abnormal red drop.

    Parameters
    ----------
    method : "ksd" or "fixed"
        "ksd"   - data-driven cutoff per group: mean + k * SD of the drop.
        "fixed" - the same absolute cutoff (e.g. 0.2 = 20 % drop) everywhere.
    k : SD multiplier for "ksd" (k = 2 is sensitive but still conservative).
    threshold : cutoff for "fixed".
    group_cols : columns defining an experimental condition; thresholds and
        summaries are computed separately for each group.

    Returns
    -------
    annotated : input rows plus metrics, thresholds and boolean flags.
    summary   : per-group counts and the fraction of cells removed.
    """
    if method not in ("ksd", "fixed"):
        raise ValueError("method must be 'ksd' or 'fixed'")
    group_cols = list(group_cols)

    d = compute_drops(df, use_percent_drop)
    thr = _thresholds(d, group_cols, method, threshold, k)
    d = d.merge(thr[group_cols + [f"threshold_{n}" for n in WINDOWS]],
                on=group_cols, how="left")

    # NaN metrics compare as False, so cells with missing readouts are kept.
    for name in WINDOWS:
        d[f"flag_{name}"] = d[f"metric_{name}"] > d[f"threshold_{name}"]
    d["flag_any"] = d[[f"flag_{n}" for n in WINDOWS]].any(axis=1)

    has_metric = d[[f"metric_{n}" for n in WINDOWS]].notna().any(axis=1)
    counts = (
        d.assign(_defined=has_metric)
        .groupby(group_cols, dropna=False)
        .agg(n_defined=("_defined", "sum"),
             n_flagged_r0_r3=("flag_r0_r3", "sum"),
             n_flagged_r3_r5=("flag_r3_r5", "sum"),
             n_flagged_any=("flag_any", "sum"))
        .reset_index()
    )
    summary = thr.merge(counts, on=group_cols, how="left")
    summary["fraction_removed"] = summary["n_flagged_any"] / summary["n_defined"].replace(0, np.nan)
    return d, summary


def filter_cells(annotated: pd.DataFrame) -> pd.DataFrame:
    """Return only the cells that passed QC."""
    return annotated.loc[~annotated["flag_any"]].copy()


def sensitivity(
    df: pd.DataFrame,
    method: str = "ksd",
    values: Iterable[float] = (1.0, 1.5, 2.0, 2.5, 3.0),
    group_cols: Sequence[str] = ("condition",),
    use_percent_drop: bool = True,
) -> pd.DataFrame:
    """Fraction of cells removed per group across a range of k (or fixed thresholds).

    Useful for justifying the chosen cutoff: a good value sits where the
    removal rate has levelled off, i.e. only clear outliers are removed.
    """
    rows = []
    for v in values:
        kwargs = {"k": v} if method == "ksd" else {"threshold": v}
        _, s = flag_cells(df, method=method, group_cols=group_cols,
                          use_percent_drop=use_percent_drop, **kwargs)
        s = s[list(group_cols) + ["fraction_removed"]].assign(parameter=v)
        rows.append(s)
    return pd.concat(rows, ignore_index=True)
