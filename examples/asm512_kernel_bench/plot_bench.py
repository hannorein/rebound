#!/usr/bin/env python3
"""Plot Newton-iteration counts vs simulation time for the kernel benchmark."""
import csv
import argparse
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CONFIGS = ["asm512", "asm512_mom", "asm512_fused", "asm512_opt"]
COLORS  = {"asm512":"#888888", "asm512_mom":"#1f77b4",
           "asm512_fused":"#2ca02c", "asm512_opt":"#d62728"}
LSTYLE  = {"asm512":":", "asm512_mom":"-",
           "asm512_fused":"-", "asm512_opt":"-"}

def load(path):
    data = defaultdict(lambda: {"t":[], "wall":[], "ctr":[], "iters":[]})
    with open(path) as f:
        for row in csv.DictReader(f):
            c = row["config"]
            data[c]["t"].append(float(row["t_yr"]))
            data[c]["wall"].append(float(row["wall_s_chunk"]))
            data[c]["ctr"].append(int(row["counter_chunk"]))
            v = row["ctr_per_step_per_lane"]
            data[c]["iters"].append(float("nan") if v == "nan" else float(v))
    return data

def plot(csv_path, png_path, setup_label):
    d = load(csv_path)
    fig, (ax_it, ax_w) = plt.subplots(1, 2, figsize=(11.5, 4.2))
    for c in CONFIGS:
        if c not in d: continue
        t = np.array(d[c]["t"])
        it = np.array(d[c]["iters"])
        w = np.array(d[c]["wall"])
        finite = ~np.isnan(it)
        ax_it.plot(t[finite], it[finite] + 1.0,
                   color=COLORS[c], linestyle=LSTYLE[c], linewidth=1.5, label=c)
        ax_w.plot(t, w * 1000.0,
                  color=COLORS[c], linestyle=LSTYLE[c], linewidth=1.5, label=c)
    ax_it.set_xlabel("simulation time [yr]")
    ax_it.set_ylabel("avg Newton iters per Kepler step (per lane)")
    ax_it.set_title(f"Newton iters: {setup_label}")
    ax_it.grid(True, alpha=0.3)
    ax_it.legend(loc="best", fontsize=9)
    all_iters = []
    for c in CONFIGS:
        if c in d:
            arr = np.array(d[c]["iters"])
            all_iters.extend((arr[~np.isnan(arr)] + 1.0).tolist())
    if all_iters:
        ymin, ymax = min(all_iters), max(all_iters)
        pad = max(0.02, 0.1 * (ymax - ymin))
        ax_it.set_ylim(ymin - pad, ymax + pad)
    ax_w.set_xlabel("simulation time [yr]")
    ax_w.set_ylabel("wall time per chunk [ms]")
    ax_w.set_title(f"Wall time: {setup_label}")
    ax_w.grid(True, alpha=0.3)
    ax_w.legend(loc="best", fontsize=9)
    fig.suptitle(f"asm512 kernel variants — {setup_label}", y=1.02)
    fig.tight_layout()
    fig.savefig(png_path, dpi=130, bbox_inches="tight")
    print(f"wrote {png_path}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--label", default="")
    a = p.parse_args()
    plot(a.csv, a.out, a.label)

if __name__ == "__main__":
    main()
