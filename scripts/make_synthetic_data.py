"""Generate a synthetic iodide-quench dataset with known dead/detached cells.

No real experimental data is distributed with this repository. The synthetic
data mimic the structure of the real pipeline output (one row per cell, red
signal at readouts r_0, r_3, r_5) and include a known set of "ejected" cells,
so the filter can be checked against ground truth.

    python scripts/make_synthetic_data.py --out data/synthetic_quench.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

CONDITIONS = ["WT_DMSO", "F508del_DMSO", "F508del_corrector"]


def make(n_per_condition: int = 300, frac_bad: float = 0.06, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for cond in CONDITIONS:
        n = n_per_condition
        r0 = rng.lognormal(mean=np.log(1000), sigma=0.35, size=n)
        drop1 = rng.normal(0.03, 0.04, n)   # healthy cells: small bleaching / exchange loss
        drop2 = rng.normal(0.05, 0.05, n)

        bad = rng.random(n) < frac_bad
        window = rng.integers(0, 2, n)      # which buffer exchange the cell is lost in
        drop1 = np.where(bad & (window == 0), rng.uniform(0.4, 0.95, n), drop1)
        drop2 = np.where(bad & (window == 1), rng.uniform(0.4, 0.95, n), drop2)

        r3 = r0 * (1 - drop1)
        r5 = r3 * (1 - drop2)
        rows.append(pd.DataFrame({
            "cell_id": [f"{cond}_{i:04d}" for i in range(n)],
            "condition": cond,
            "well": rng.choice(["B02", "B03", "B04", "C02", "C03", "C04"], n),
            "r_0": r0.round(2), "r_3": r3.round(2), "r_5": r5.round(2),
            "true_bad": bad,
        }))
    return pd.concat(rows, ignore_index=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("data/synthetic_quench.csv"))
    p.add_argument("-n", type=int, default=300, help="cells per condition")
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    make(a.n, seed=a.seed).to_csv(a.out, index=False)
    print(f"Wrote {a.out}")
