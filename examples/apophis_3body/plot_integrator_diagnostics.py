#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def read_diag_csv(path):
    t = []
    dt = []
    rel_e = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        required = {"t", "dt", "relE"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"{path} missing required columns: {sorted(required)}")
        for row in reader:
            t.append(float(row["t"]))
            dt.append(float(row["dt"]))
            rel_e.append(abs(float(row["relE"])))
    return t, dt, rel_e


def main():
    parser = argparse.ArgumentParser(
        description="Overlay WHFAST and IAS15 diagnostics from REBOUND CSV outputs."
    )
    parser.add_argument("--whfast", default="whfast.csv", help="Path to WHFAST CSV")
    parser.add_argument("--ias15", default="IAS15.csv", help="Path to IAS15 CSV")
    parser.add_argument("--out", default="integrator_compare.png", help="Output figure file")
    parser.add_argument(
        "--time-years",
        action="store_true",
        help="Convert time axis from days to years (365.25 d/yr).",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    whfast_path = Path(args.whfast)
    ias15_path = Path(args.ias15)

    # Resolve relative paths from the script directory so execution works from any cwd.
    if not whfast_path.is_absolute():
        whfast_path = script_dir / whfast_path
    if not ias15_path.is_absolute():
        ias15_path = script_dir / ias15_path

    t_w, dt_w, e_w = read_diag_csv(whfast_path)
    t_i, dt_i, e_i = read_diag_csv(ias15_path)

    if args.time_years:
        scale = 1.0 / 365.25
        t_w = [x * scale for x in t_w]
        t_i = [x * scale for x in t_i]
        xlabel = "time [yr]"
    else:
        xlabel = "time [day]"

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    ax1.plot(t_w, dt_w, label="WHFAST", lw=1.5)
    ax1.plot(t_i, dt_i, label="IAS15", lw=1.5)
    ax1.set_ylabel("dt [day]")
    ax1.set_title("Timestep Size vs Time")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    ax2.plot(t_w, e_w, label="WHFAST", lw=1.5)
    ax2.plot(t_i, e_i, label="IAS15", lw=1.5)
    ax2.set_yscale("log")
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel("|relE|")
    ax2.set_title("Global Relative Energy Error vs Time")
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend()

    fig.tight_layout()
    # fig.savefig(args.out, dpi=150)
    plt.show(block=True)
    # print(f"Saved plot to {args.out}")


if __name__ == "__main__":
    main()
