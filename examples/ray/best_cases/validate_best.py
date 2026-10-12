#!/usr/bin/env python3
"""Validate the retained cases with a fixed-duration timestep scan."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from search_fixed_hierarchy import DEFAULT_OUTPUT, evaluate, kepler_period

import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--divisors", type=int, nargs="+", default=(10, 20, 40, 80, 160))
    parser.add_argument("--samples", type=int, default=1001)
    parser.add_argument("--duration-multiplier", type=float, default=1.0)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.duration_multiplier <= 0.0:
        raise ValueError("duration multiplier must be positive")

    rank_dirs = sorted(path for path in args.input.glob("rank_*"))
    all_rows = []
    for rank_dir in rank_dirs:
        with (rank_dir / "parameters.json").open() as stream:
            base = json.load(stream)
        base["duration"] *= args.duration_multiplier
        outer_period = kepler_period(base["a_outer"], sum(base["masses"]))
        rows = []
        for divisor in args.divisors:
            config = dict(base)
            config["dt_divisor"] = divisor
            config["dt"] = config["reference_period"] / divisor
            config["n_steps"] = int(math.ceil(config["duration"] / config["dt"]))
            result = evaluate(config, args.samples)
            row = {
                "rank": rank_dir.name,
                "id": base["id"],
                "dt_divisor": divisor,
                "dt_over_reference_period": 1.0 / divisor,
                "n_steps": config["n_steps"],
                "duration": config["duration"],
                "duration_over_outer_period": config["duration"] / outer_period,
                "duration_multiplier": args.duration_multiplier,
                "whfast_max_error": result["whfast_max_error"],
                "whfast_hj_max_error": result["whfast_hj_max_error"],
                "improvement": result["improvement"],
            }
            rows.append(row)
            all_rows.append(row)
            print(
                f"{rank_dir.name} P/{divisor:<3d}: "
                f"WH={row['whfast_max_error']:.3e} "
                f"HJ={row['whfast_hj_max_error']:.3e} "
                f"ratio={row['improvement']:.3e}"
            )

        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3))
        x = np.asarray([row["dt_over_reference_period"] for row in rows])
        wh = np.asarray([row["whfast_max_error"] for row in rows])
        hj = np.asarray([row["whfast_hj_max_error"] for row in rows])
        ratio = wh / hj
        order = np.argsort(x)
        axes[0].loglog(x[order], wh[order], "o-", label="WHFast", color="#d95f02")
        axes[0].loglog(x[order], hj[order], "s-", label="fixed HJ", color="#1b9e77")
        axes[0].set_xlabel(r"$\Delta t/P_{\rm ref}$")
        axes[0].set_ylabel("maximum relative energy error")
        axes[0].legend()
        axes[1].loglog(x[order], ratio[order], "o-", color="#3b5bdb")
        axes[1].axhline(1.0e4, color="black", ls="--", lw=1, label="10,000x target")
        axes[1].set_xlabel(r"$\Delta t/P_{\rm ref}$")
        axes[1].set_ylabel("WHFast error / fixed-HJ error")
        axes[1].legend()
        for ax in axes:
            ax.grid(True, which="both", alpha=0.25)
        fig.suptitle(f"Timestep validation: {rank_dir.name} ({base['id']})")
        fig.tight_layout()
        fig.savefig(rank_dir / "timestep_validation.png", dpi=180)
        plt.close(fig)

    if all_rows:
        with (args.input / "timestep_validation.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(all_rows[0]))
            writer.writeheader()
            writer.writerows(all_rows)


if __name__ == "__main__":
    main()
