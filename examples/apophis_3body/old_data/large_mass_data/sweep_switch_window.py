#!/usr/bin/env python3
"""Sweep symmetric switch windows around closest approach for the Apophis example."""

# python3 examples/apophis_3body/sweep_switch_window.py \
#   --summary-csv examples/apophis_3body/my_sweep.csv \
#   --plot examples/apophis_3body/my_sweep.png
from __future__ import annotations

import argparse
from bisect import bisect_left
import csv
import math
import sys
import tempfile
import warnings
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[3]
sys.path.insert(0, str(REPO_ROOT))

import examples.apophis_3body.old_data.large_mass_data.run_saved_snapshot_compare as compare

SWEEP_MODES = (
    {
        "mode": "whfast_fake_hj",
        "label": "WHFast -> fake HJ -> WHFast",
        "runner": compare.run_fake_hj_window_case,
    },
    {
        "mode": "ias15_whfast_fake_hj",
        "label": "IAS15 -> WHFast+fake HJ -> IAS15",
        "runner": compare.run_ias15_whfast_fake_hj_case,
    },
)


def earth_asteroid_distance(sim: compare.rebound.Simulation) -> float:
    """Return the Earth-Asteroid separation in AU."""
    particles = sim.particles
    dx = particles["asteroid"].x - particles["earth"].x
    dy = particles["asteroid"].y - particles["earth"].y
    dz = particles["asteroid"].z - particles["earth"].z
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def load_base_sim(snapshot_path: Path) -> compare.rebound.Simulation:
    """Load the saved snapshot and restore Python-side function pointers."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="You have to reset function pointers after creating a reb_simulation struct with a binary file.",
            category=RuntimeWarning,
        )
        base_sim = compare.rebound.Simulation(str(snapshot_path))
    compare.rebound.clibrebound.reb_simulation_reset_function_pointers(compare.byref(base_sim))
    for particle, key in zip(base_sim.particles[:3], compare.STANDARD_ORDER):
        particle.hash = key
    return base_sim


def find_closest_approach(
    base_sim: compare.rebound.Simulation,
    duration_days: float,
    whfast_dt: float,
) -> tuple[float, float]:
    """Return the time and distance of closest Earth-Asteroid approach."""
    sim = compare.prepare_sim(base_sim, "ias15", whfast_dt)
    t_stop = sim.t + duration_days
    min_distance = earth_asteroid_distance(sim)
    min_time = sim.t

    while sim.t < t_stop:
        sim.step()
        distance = earth_asteroid_distance(sim)
        if distance < min_distance:
            min_distance = distance
            min_time = sim.t

    return min_time, min_distance


def distances_at_times(
    base_sim: compare.rebound.Simulation,
    times: list[float],
    whfast_dt: float,
) -> dict[float, float]:
    """Evaluate the Earth-Asteroid distance at exact target times with IAS15."""
    if not times:
        return {}

    sim = compare.prepare_sim(base_sim, "ias15", whfast_dt)
    distances: dict[float, float] = {}

    for target_time in sorted(set(times)):
        sim.integrate(target_time, exact_finish_time=1)
        distances[target_time] = earth_asteroid_distance(sim)

    return distances


def linspace(start: float, stop: float, num: int) -> list[float]:
    """Small numpy-free linspace helper."""
    if num <= 0:
        raise ValueError("Number of sweep points must be positive.")
    if num == 1:
        return [start]
    step = (stop - start) / (num - 1)
    return [start + i * step for i in range(num)]


def write_summary_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    """Write the sweep summary table."""
    fieldnames = [
        "mode",
        "case",
        "desired_half_window_days",
        "desired_switch_start_days",
        "desired_switch_end_days",
        "closest_approach_days",
        "closest_approach_distance_AU",
        "actual_switch_start_days",
        "actual_switch_end_days",
        "actual_start_offset_before_ca_days",
        "actual_end_offset_after_ca_days",
        "desired_switch_start_distance_AU",
        "desired_switch_end_distance_AU",
        "actual_switch_start_distance_AU",
        "actual_switch_end_distance_AU",
        "steps_taken",
        "final_relE",
        "max_relE",
    ]

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def interpolate_reference_distance(offset_days: float, offsets: list[float], distances: list[float]) -> float:
    """Linearly interpolate the reference inbound distance for a time-before-CA offset."""
    if not offsets:
        raise ValueError("Reference offset grid is empty.")
    if len(offsets) == 1:
        return distances[0]
    if offset_days <= offsets[0]:
        return distances[0]
    if offset_days >= offsets[-1]:
        return distances[-1]

    upper = bisect_left(offsets, offset_days)
    lower = upper - 1
    x0 = offsets[lower]
    x1 = offsets[upper]
    y0 = distances[lower]
    y1 = distances[upper]
    if x1 == x0:
        return y1
    fraction = (offset_days - x0) / (x1 - x0)
    return y0 + fraction * (y1 - y0)


def make_plot(path: Path, rows: list[dict[str, float | str]]) -> None:
    """Plot error and total steps against switch time with a distance reference axis."""
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "matplotlib is required to write the sweep plot. "
            "Re-run with --skip-plot to generate only the summary CSV."
        ) from exc

    valid_rows = [
        row for row in rows if row["actual_switch_start_days"] != "" and row["actual_switch_start_distance_AU"] != ""
    ]
    if not valid_rows:
        raise RuntimeError("No valid sweep rows to plot.")

    reference_map = {
        float(row["actual_start_offset_before_ca_days"]): float(row["actual_switch_start_distance_AU"])
        for row in valid_rows
    }
    reference_offsets = sorted(reference_map)
    reference_distances = [reference_map[offset] for offset in reference_offsets]

    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    ax_steps = ax.twinx()

    style_map = {
        "whfast_fake_hj": {"color": "#1f77b4", "marker": "o"},
        "ias15_whfast_fake_hj": {"color": "#d62728", "marker": "s"},
    }
    label_map = {entry["mode"]: entry["label"] for entry in SWEEP_MODES}

    for mode in label_map:
        mode_rows = [row for row in valid_rows if row["mode"] == mode]
        if not mode_rows:
            continue
        rows_by_time = sorted(mode_rows, key=lambda row: float(row["actual_start_offset_before_ca_days"]))
        style = style_map.get(mode, {"marker": "o"})

        ax.plot(
            [float(row["actual_start_offset_before_ca_days"]) for row in rows_by_time],
            [float(row["max_relE"]) for row in rows_by_time],
            label=label_map[mode],
            lw=1.5,
            **style,
        )
        ax_steps.plot(
            [float(row["actual_start_offset_before_ca_days"]) for row in rows_by_time],
            [float(row["steps_taken"]) for row in rows_by_time],
            lw=1.2,
            ls="--",
            alpha=0.8,
            color=style.get("color"),
        )

    ax.set_yscale("log")
    ax.set_xlabel("actual switch start before closest approach [day]")
    ax.set_ylabel("max |relE|")
    ax.set_title("Maximum Relative Energy Error and Total Steps vs Switch Time")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    ax_steps.set_ylabel("total steps (dashed)")

    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    bottom_ticks = [tick for tick in ax.get_xticks() if reference_offsets[0] <= tick <= reference_offsets[-1]]
    if bottom_ticks:
        ax_top.set_xticks(bottom_ticks)
        ax_top.set_xticklabels(
            [
                f"{interpolate_reference_distance(tick, reference_offsets, reference_distances):.4g}"
                for tick in bottom_ticks
            ]
        )
    ax_top.set_xlabel("IAS15 reference inbound switch distance [AU]")

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def build_row(
    mode: str,
    case_index: int,
    closest_approach_time: float,
    closest_approach_distance: float,
    half_window_days: float,
    summary: dict[str, float],
    distance_lookup: dict[float, float],
) -> dict[str, float | str]:
    """Combine run diagnostics with reference distances for one sweep point."""
    desired_start = closest_approach_time - half_window_days
    desired_end = closest_approach_time + half_window_days
    actual_start = summary.get("switch_start_actual")
    actual_end = summary.get("switch_end_actual")

    row: dict[str, float | str] = {
        "mode": mode,
        "case": case_index,
        "desired_half_window_days": half_window_days,
        "desired_switch_start_days": desired_start,
        "desired_switch_end_days": desired_end,
        "closest_approach_days": closest_approach_time,
        "closest_approach_distance_AU": closest_approach_distance,
        "actual_switch_start_days": "" if actual_start is None else actual_start,
        "actual_switch_end_days": "" if actual_end is None else actual_end,
        "actual_start_offset_before_ca_days": "" if actual_start is None else closest_approach_time - actual_start,
        "actual_end_offset_after_ca_days": "" if actual_end is None else actual_end - closest_approach_time,
        "desired_switch_start_distance_AU": distance_lookup[desired_start],
        "desired_switch_end_distance_AU": distance_lookup[desired_end],
        "actual_switch_start_distance_AU": "" if actual_start is None else distance_lookup[actual_start],
        "actual_switch_end_distance_AU": "" if actual_end is None else distance_lookup[actual_end],
        "steps_taken": int(summary["steps_taken"]),
        "final_relE": summary["final_relE"],
        "max_relE": summary["max_relE"],
    }
    return row


def print_summary(rows: list[dict[str, float | str]]) -> None:
    """Print a compact text table to stdout."""
    headers = (
        "mode",
        "case",
        "half_win[d]",
        "start_before_CA[d]",
        "start_dist[AU]",
        "max|relE|",
        "final|relE|",
    )
    print(
        f"{headers[0]:>22}  {headers[1]:>4}  {headers[2]:>11}  {headers[3]:>17}  "
        f"{headers[4]:>14}  {headers[5]:>12}  {headers[6]:>12}"
    )
    print(
        f"{'-' * 22}  {'-' * 4}  {'-' * 11}  {'-' * 17}  "
        f"{'-' * 14}  {'-' * 12}  {'-' * 12}"
    )

    for row in rows:
        start_before_ca = row["actual_start_offset_before_ca_days"]
        start_distance = row["actual_switch_start_distance_AU"]
        start_before_ca_text = "-" if start_before_ca == "" else f"{float(start_before_ca):.6g}"
        start_distance_text = "-" if start_distance == "" else f"{float(start_distance):.6g}"
        print(
            f"{str(row['mode']):>22}  "
            f"{int(row['case']):>4d}  "
            f"{float(row['desired_half_window_days']):>11.6g}  "
            f"{start_before_ca_text:>17}  "
            f"{start_distance_text:>14}  "
            f"{float(row['max_relE']):>12.3e}  "
            f"{float(row['final_relE']):>12.3e}"
        )


def main() -> None:
    default_snapshot = SCRIPT_DIR.parents[1] / "find_close_encounter" / "apophis_1yr_before.bin"

    parser = argparse.ArgumentParser(
        description=(
            "Sweep symmetric switch windows around the closest Earth-Asteroid approach and "
            "overlay WHFast -> fake-HJ -> WHFast with IAS15 -> WHFast+fake-HJ -> IAS15."
        )
    )
    parser.add_argument("--snapshot", type=Path, default=default_snapshot, help="Saved REBOUND binary snapshot.")
    parser.add_argument(
        "--duration-days",
        type=float,
        default=365.25 * 2.0,
        help="Integration span after the saved snapshot.",
    )
    parser.add_argument(
        "--whfast-dt",
        type=float,
        default=10, 
        help="WHFast timestep in days inside the switched window.",
    )
    parser.add_argument(
        "--num-windows",
        type=int,
        default=15,
        help="Number of symmetric switch windows to test.",
    )
    parser.add_argument(
        "--min-half-window-days",
        type=float,
        default=None,
        help="Smallest half-width around closest approach to test. Defaults to WHFast dt.",
    )
    parser.add_argument(
        "--max-half-window-days",
        type=float,
        default=None,
        help="Largest half-width around closest approach to test. Defaults to min(180 d, available window).",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=compare.SCRIPT_DIR / "switch_window_sweep_summary.csv",
        help="Output CSV containing one row per experiment run.",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=compare.SCRIPT_DIR / "switch_window_sweep.png",
        help="Output plot file.",
    )
    parser.add_argument(
        "--skip-plot",
        action="store_true",
        help="Skip figure generation and write only the summary CSV.",
    )
    args = parser.parse_args()

    snapshot_path = args.snapshot.resolve()
    if not snapshot_path.exists():
        raise FileNotFoundError(f"Snapshot file not found: {snapshot_path}")
    if args.num_windows <= 0:
        raise ValueError("--num-windows must be positive.")

    base_sim = load_base_sim(snapshot_path)
    t_start = base_sim.t
    t_stop = t_start + args.duration_days

    closest_approach_time, closest_approach_distance = find_closest_approach(
        base_sim,
        args.duration_days,
        args.whfast_dt,
    )
    available_before = closest_approach_time - t_start
    available_after = t_stop - closest_approach_time
    available_half_window = min(available_before, available_after)

    min_half_window = args.min_half_window_days
    if min_half_window is None:
        min_half_window = args.whfast_dt
    max_half_window = args.max_half_window_days
    if max_half_window is None:
        max_half_window = min(180.0, 0.95 * available_half_window)

    if min_half_window <= 0.0:
        raise ValueError("Minimum half window must be positive.")
    if max_half_window <= 0.0:
        raise ValueError("Maximum half window must be positive.")
    if min_half_window > max_half_window:
        raise ValueError("Minimum half window must be <= maximum half window.")
    if max_half_window >= available_half_window:
        raise ValueError(
            f"Maximum half window {max_half_window} d exceeds the available symmetric interval "
            f"{available_half_window} d around closest approach."
        )

    half_windows = linspace(min_half_window, max_half_window, args.num_windows)

    summaries: list[tuple[str, int, float, dict[str, float]]] = []
    print(
        f"Closest approach at t = {closest_approach_time:.10f} d "
        f"with d = {closest_approach_distance:.12g} AU"
    )
    print(
        f"Sweeping {len(half_windows)} symmetric windows from {half_windows[0]:.6g} d "
        f"to {half_windows[-1]:.6g} d around closest approach."
    )

    with tempfile.TemporaryDirectory(prefix="apophis_switch_sweep_") as tmpdir:
        tmpdir_path = Path(tmpdir)
        for mode_entry in SWEEP_MODES:
            print(f"Running {mode_entry['label']}...", flush=True)
            for case_index, half_window_days in enumerate(half_windows, start=1):
                switch_start_yr = (closest_approach_time - half_window_days - t_start) / compare.YEAR_DAYS
                switch_end_yr = (closest_approach_time + half_window_days - t_start) / compare.YEAR_DAYS
                output_csv = tmpdir_path / f"{mode_entry['mode']}_case_{case_index:03d}.csv"
                summary = mode_entry["runner"](
                    base_sim,
                    output_csv,
                    args.duration_days,
                    args.whfast_dt,
                    switch_start_yr,
                    switch_end_yr,
                )
                summaries.append((mode_entry["mode"], case_index, half_window_days, summary))

    query_times = [closest_approach_time]
    for _, _, half_window_days, summary in summaries:
        desired_start = closest_approach_time - half_window_days
        desired_end = closest_approach_time + half_window_days
        query_times.extend([desired_start, desired_end])
        actual_start = summary.get("switch_start_actual")
        actual_end = summary.get("switch_end_actual")
        if actual_start is not None:
            query_times.append(actual_start)
        if actual_end is not None:
            query_times.append(actual_end)

    distance_lookup = distances_at_times(base_sim, query_times, args.whfast_dt)

    rows = [
        build_row(
            mode,
            case_index,
            closest_approach_time,
            closest_approach_distance,
            half_window_days,
            summary,
            distance_lookup,
        )
        for mode, case_index, half_window_days, summary in summaries
    ]

    write_summary_csv(args.summary_csv.resolve(), rows)
    print_summary(rows)
    print(f"Saved summary CSV to {args.summary_csv.resolve()}")
    if args.skip_plot:
        print("Skipped plot generation (--skip-plot).")
    else:
        make_plot(args.plot.resolve(), rows)
        print(f"Saved plot to {args.plot.resolve()}")


if __name__ == "__main__":
    main()
