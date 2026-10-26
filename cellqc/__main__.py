"""Command-line interface: python -m cellqc INPUT [options]"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .filters import filter_cells, flag_cells, sensitivity


def read_table(path: Path, sheet: str | None) -> pd.DataFrame:
    if path.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(path, sheet_name=sheet or 0)
    return pd.read_csv(path)


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(
        prog="cellqc",
        description="Flag and remove dead/detached cells from iodide-quench single-cell data.",
    )
    p.add_argument("input", type=Path, help="CSV or Excel file with columns r_0, r_3, r_5 and a condition column")
    p.add_argument("--sheet", default=None, help="Excel sheet name (default: first sheet)")
    p.add_argument("--method", choices=["ksd", "fixed"], default="ksd")
    p.add_argument("-k", type=float, default=2.0, help="SD multiplier for --method ksd (default 2)")
    p.add_argument("--threshold", type=float, default=0.2, help="cutoff for --method fixed (default 0.2)")
    p.add_argument("--group", nargs="+", default=["condition"], help="grouping column(s)")
    p.add_argument("--absolute", action="store_true", help="use absolute drop instead of fractional drop")
    p.add_argument("--sensitivity", action="store_true", help="also write a sensitivity sweep and plot")
    p.add_argument("-o", "--outdir", type=Path, default=Path("results"))
    a = p.parse_args(argv)

    df = read_table(a.input, a.sheet)
    annotated, summary = flag_cells(
        df, method=a.method, k=a.k, threshold=a.threshold,
        group_cols=a.group, use_percent_drop=not a.absolute,
    )
    a.outdir.mkdir(parents=True, exist_ok=True)
    stem = a.input.stem
    annotated.to_csv(a.outdir / f"{stem}_annotated.csv", index=False)
    filter_cells(annotated).to_csv(a.outdir / f"{stem}_filtered.csv", index=False)
    summary.to_csv(a.outdir / f"{stem}_summary.csv", index=False)

    n, bad = len(annotated), int(annotated["flag_any"].sum())
    print(f"{bad}/{n} cells flagged ({bad / n:.1%}) using method={a.method}")
    print(summary[a.group + ["threshold_r0_r3", "threshold_r3_r5", "fraction_removed"]].to_string(index=False))

    if a.sensitivity:
        values = [1.0, 1.5, 2.0, 2.5, 3.0] if a.method == "ksd" else [0.1, 0.15, 0.2, 0.25, 0.3]
        sweep = sensitivity(df, method=a.method, values=values,
                            group_cols=a.group, use_percent_drop=not a.absolute)
        sweep.to_csv(a.outdir / f"{stem}_sensitivity.csv", index=False)

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(5, 3.5))
        for key, g in sweep.groupby(a.group):
            label = ", ".join(map(str, key)) if isinstance(key, tuple) else str(key)
            ax.plot(g["parameter"], g["fraction_removed"], marker="o", label=label)
        ax.set_xlabel("k (SD multiplier)" if a.method == "ksd" else "fixed threshold")
        ax.set_ylabel("fraction of cells removed")
        ax.legend(title=", ".join(a.group), fontsize=8)
        fig.tight_layout()
        fig.savefig(a.outdir / f"{stem}_sensitivity.png", dpi=200)
        print("Wrote sensitivity sweep and plot")

    print(f"Results in {a.outdir.resolve()}")


if __name__ == "__main__":
    main()
