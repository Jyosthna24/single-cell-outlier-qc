# single-cell-outlier-qc

![tests](https://github.com/Jyosthna24/single-cell-outlier-qc/actions/workflows/tests.yml/badge.svg)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)

**Statistical quality control that removes dead and detached cells from high-content imaging data before downstream analysis.**

I wrote this for my PhD (UCL, Pharmacology & Physiology). There it runs as a QC step in a single-cell image-analysis pipeline for **CFTR drug discovery**. Each experiment produces thousands of cells per plate. A few percent of them die or are washed off during buffer exchange. Their signal collapses for reasons that have nothing to do with the biology being measured, and if they stay in the data they bias every downstream estimate.

---

## The problem

In the iodide-quench assay, cells express a red reference fluorophore (mCherry) whose signal should stay roughly **constant** over the experiment. A large, sudden drop in red signal means the cell has been lost, not that the drug did something.

| | Healthy cell | Dead / detached cell |
|---|---|---|
| Red signal r₀ → r₃ → r₅ | small, gradual loss (bleaching) | sudden 40–95 % drop |
| Keep in analysis? | ✅ | ❌ |

Removing these cells by eye is slow and impossible to reproduce. This package replaces that step with an objective, documented rule.

## Method

For every cell, I compute the **fractional drop** in red signal over two windows. A fractional drop is scale-free, so it can be compared across cells with very different expression levels:

$$
D_{0\to3} = \frac{r_0 - r_3}{r_0}, \qquad D_{3\to5} = \frac{r_3 - r_5}{r_3}
$$

A cell is removed if **either** drop exceeds a threshold *T*:

$$
\text{remove} \iff D_{0\to3} > T_{0\to3} \;\lor\; D_{3\to5} > T_{3\to5}
$$

There are two ways to set *T*:

| Method | Threshold | When to use |
|---|---|---|
| `ksd` (data-driven) | $T = \mu + k\sigma$ of the drop, computed **separately for each experimental condition** | When conditions differ in baseline behaviour. k = 2 is sensitive but conservative. |
| `fixed` | $T = 0.2$ (a 20 % drop) everywhere | When you want the same absolute rule applied to every plate and condition |

Edge cases are handled explicitly. A zero baseline gives `NaN`, not `inf`, and a cell with a missing readout is kept rather than silently dropped.

### Choosing the cutoff: sensitivity analysis

A threshold should be justified, not picked. Run with `--sensitivity` to sweep *k* (or *T*) and plot the fraction of cells removed in each condition. A good cutoff sits where that curve flattens out, the point where only clear outliers are being removed.

### Validation on ground truth

`scripts/make_synthetic_data.py` simulates a plate with a **known** set of dead cells. The test suite checks that the filter recovers them with **recall > 95 % and precision > 90 %**.

## Usage

```bash
pip install -r requirements.txt

# 1. generate example data (no real experimental data is shipped)
python scripts/make_synthetic_data.py --out data/synthetic_quench.csv

# 2. data-driven filter, per condition, with sensitivity sweep
python -m cellqc data/synthetic_quench.csv --method ksd -k 2 --sensitivity

# or a fixed 20 % cutoff
python -m cellqc data/synthetic_quench.csv --method fixed --threshold 0.2

# real pipeline output (Excel) works too
python -m cellqc full_quench_results.xlsx --sheet full_results
```

Input: one row per cell, with red-signal columns `r_0`, `r_3`, `r_5` and a grouping column (default `condition`; use `--group` to change it).

Outputs, written to `results/`:

| File | Contents |
|---|---|
| `*_annotated.csv` | every cell, with its drop metrics, thresholds and flags |
| `*_filtered.csv` | only the cells that passed QC |
| `*_summary.csv` | per condition: n, mean, SD, threshold, cells removed, fraction removed |
| `*_sensitivity.csv/.png` | fraction removed across the parameter sweep |

As a library:

```python
from cellqc import flag_cells, filter_cells
annotated, summary = flag_cells(df, method="ksd", k=2, group_cols=["condition"])
clean = filter_cells(annotated)
```

## Repository layout

```
cellqc/            core library (filters.py) and command-line interface (__main__.py)
scripts/           synthetic data generator with ground-truth labels
tests/             pytest suite, run automatically by GitHub Actions on every push
```

## Skills demonstrated
Outlier detection · per-group robust thresholds · sensitivity analysis · validation against simulated ground truth · reproducible CLI tooling · unit testing and CI · pandas and NumPy

---

*Part of my PhD work on high-content screening for CFTR modulators. See also my [GitHub profile](https://github.com/Jyosthna24) and [portfolio](https://jyosthna24.github.io).*
