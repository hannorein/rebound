#!/usr/bin/env python3
"""Compare WHFast and IAS15 from the saved 1-year-before-encounter snapshot."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import sys
import sysconfig
from ctypes import byref

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(REPO_ROOT))

# Keep matplotlib cache writes inside the repo workspace for portability.
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".mplconfig"))

suffix = sysconfig.get_config_var("EXT_SUFFIX") or ".so"
expected_lib = REPO_ROOT / f"librebound{suffix}"
source_lib = REPO_ROOT / "src" / "librebound.so"
if not expected_lib.exists() and source_lib.exists():
    expected_lib.symlink_to(source_lib.relative_to(REPO_ROOT))

import matplotlib.pyplot as plt
import rebound


GAUSSIAN_K = 0.01720209895
YEAR_DAYS = 365.25
STANDARD_ORDER = ("sun", "earth", "asteroid")
FAKE_HJ_ORDER = ("earth", "asteroid", "sun")
STANDARD_ORDER_NAME = "sun-earth-asteroid"
FAKE_HJ_ORDER_NAME = "earth-asteroid-sun"


def clone_sim_with_order(
    source_sim: rebound.Simulation,
    order: tuple[str, ...],
    integrator: str,
    whfast_dt: float,
) -> rebound.Simulation:
    """Rebuild a simulation from the current inertial state with a new particle order."""
    new_sim = rebound.Simulation()
    new_sim.t = source_sim.t
    new_sim.softening = source_sim.softening
    new_sim.testparticle_type = source_sim.testparticle_type
    for key in order:
        new_sim.add(source_sim.particles[key].copy())
        new_sim.particles[-1].hash = key
    configure_sim(new_sim, integrator, whfast_dt)
    return new_sim


def prepare_sim(
    base_sim: rebound.Simulation,
    integrator: str,
    whfast_dt: float,
    order: tuple[str, ...] = STANDARD_ORDER,
) -> rebound.Simulation:
    """Create a simulation branch for the requested experiment."""
    if integrator == "ias15" and order == STANDARD_ORDER:
        sim = base_sim.copy()
        configure_sim(sim, integrator, whfast_dt)
        return sim
    return clone_sim_with_order(base_sim, order, integrator, whfast_dt)


def configure_sim(sim: rebound.Simulation, integrator: str, whfast_dt: float) -> None:
    """Match the existing C example settings for each integrator."""
    sim.N_active = 3
    sim.force_is_velocity_dependent = 0
    sim.G = GAUSSIAN_K * GAUSSIAN_K

    if integrator == "ias15":
        # The snapshot was saved from an IAS15 run, so keep its state if possible.
        if sim.integrator != "ias15":
            sim.reset_integrator()
            sim.integrator = "ias15"
        return

    if integrator == "whfast":
        sim.reset_integrator()
        sim.integrator = "whfast"
        sim.dt = whfast_dt
        sim.ri_whfast.safe_mode = 1
        sim.ri_whfast.corrector = 11
        sim.ri_whfast.coordinates = "jacobi"
        return

    raise ValueError(f"Unsupported integrator {integrator!r}")


def run_case(
    base_sim: rebound.Simulation,
    integrator: str,
    output_csv: Path,
    duration_days: float,
    whfast_dt: float,
) -> dict[str, float]:
    sim = prepare_sim(base_sim, integrator, whfast_dt)

    if sim.N < 3:
        raise RuntimeError(f"Snapshot must contain at least 3 particles; found N={sim.N}")

    t_start = sim.t
    t_stop = t_start + duration_days
    e_initial = sim.energy()
    step_count = 0
    max_abs_rel_e = 0.0

    with output_csv.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["step", "t", "t_rel", "dt", "relE", "ordering"])
        while sim.t < t_stop:
            remaining = t_stop - sim.t
            if remaining <= 0.0:
                break

            # Keep WHFast on a fixed step grid. For IAS15, allow the last
            # integration interval to shrink to the requested stop time.
            if integrator != "whfast" and sim.dt > remaining:
                sim.dt = remaining

            sim.step()
            step_count += 1
            rel_e = (sim.energy() - e_initial) / e_initial
            max_abs_rel_e = max(max_abs_rel_e, abs(rel_e))
            writer.writerow(
                [
                    step_count,
                    f"{sim.t:.16e}",
                    f"{(sim.t - t_start):.16e}",
                    f"{sim.dt_last_done:.16e}",
                    f"{rel_e:.16e}",
                    STANDARD_ORDER_NAME,
                ]
            )

    e_final = sim.energy()
    return {
        "t_start": t_start,
        "t_stop": sim.t,
        "steps_taken": float(step_count),
        "final_relE": abs((e_final - e_initial) / e_initial),
        "max_relE": max_abs_rel_e,
    }


def run_fake_hj_window_case(
    base_sim: rebound.Simulation,
    output_csv: Path,
    duration_days: float,
    whfast_dt: float,
    switch_start_yr: float,
    switch_end_yr: float,
) -> dict[str, float]:
    t_start = base_sim.t
    t_stop = t_start + duration_days
    switch_start = t_start + switch_start_yr * YEAR_DAYS
    switch_end = t_start + switch_end_yr * YEAR_DAYS

    if not (t_start <= switch_start < switch_end <= t_stop):
        raise ValueError(
            "Switch window must satisfy start <= switch_start < switch_end <= stop. "
            f"Got start={switch_start_yr} yr, end={switch_end_yr} yr for duration={duration_days / YEAR_DAYS} yr."
        )

    sim = prepare_sim(base_sim, "whfast", whfast_dt)
    e_initial = sim.energy()
    step_count = 0
    max_abs_rel_e = 0.0
    current_order = STANDARD_ORDER
    current_order_name = STANDARD_ORDER_NAME
    actual_switch_start = None
    actual_switch_end = None

    with output_csv.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["step", "t", "t_rel", "dt", "relE", "ordering"])

        while sim.t < t_stop:
            desired_order = STANDARD_ORDER
            desired_order_name = STANDARD_ORDER_NAME
            if switch_start <= sim.t < switch_end:
                desired_order = FAKE_HJ_ORDER
                desired_order_name = FAKE_HJ_ORDER_NAME

            if desired_order != current_order:
                sim = prepare_sim(sim, "whfast", whfast_dt, order=desired_order)
                current_order = desired_order
                current_order_name = desired_order_name
                if desired_order == FAKE_HJ_ORDER and actual_switch_start is None:
                    actual_switch_start = sim.t
                if desired_order == STANDARD_ORDER and actual_switch_end is None and sim.t >= switch_end:
                    actual_switch_end = sim.t

            sim.dt = whfast_dt
            sim.step()
            step_count += 1
            rel_e = (sim.energy() - e_initial) / e_initial
            max_abs_rel_e = max(max_abs_rel_e, abs(rel_e))
            writer.writerow(
                [
                    step_count,
                    f"{sim.t:.16e}",
                    f"{(sim.t - t_start):.16e}",
                    f"{sim.dt_last_done:.16e}",
                    f"{rel_e:.16e}",
                    current_order_name,
                ]
            )

    e_final = sim.energy()
    return {
        "t_start": t_start,
        "t_stop": sim.t,
        "steps_taken": float(step_count),
        "final_relE": abs((e_final - e_initial) / e_initial),
        "max_relE": max_abs_rel_e,
        "switch_start_actual": actual_switch_start,
        "switch_end_actual": actual_switch_end,
    }


def switch_time_for_plot(summary: dict[str, float], key: str, time_years: bool) -> float:
    """Convert an absolute switch time in the summary to plot coordinates."""
    offset = summary[key] - summary["t_start"]
    if time_years:
        return offset / YEAR_DAYS
    return offset


def read_diag_csv(path: Path):
    t_rel = []
    dt = []
    rel_e = []
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        required = {"t_rel", "dt", "relE"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"{path} missing required columns: {sorted(required)}")
        for row in reader:
            t_rel.append(float(row["t_rel"]))
            dt.append(float(row["dt"]))
            rel_e.append(abs(float(row["relE"])))
    return t_rel, dt, rel_e


def make_plot(
    whfast_csv: Path,
    ias15_csv: Path,
    fake_hj_csv: Path,
    output_png: Path,
    time_years: bool,
    shade_start: float,
    shade_end: float,
) -> None:
    t_w, dt_w, e_w = read_diag_csv(whfast_csv)
    t_i, dt_i, e_i = read_diag_csv(ias15_csv)
    t_h, dt_h, e_h = read_diag_csv(fake_hj_csv)

    if time_years:
        scale = 1.0 / YEAR_DAYS
        t_w = [x * scale for x in t_w]
        t_i = [x * scale for x in t_i]
        t_h = [x * scale for x in t_h]
        xlabel = "time since restart [yr]"
    else:
        xlabel = "time since restart [day]"

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    ax1.axvspan(shade_start, shade_end, color="0.85", alpha=0.5)
    ax1.plot(t_w, dt_w, label="WHFast", lw=1.5)
    ax1.plot(t_i, dt_i, label="IAS15", lw=1.5)
    ax1.plot(t_h, dt_h, label="WHFast fake-HJ", lw=1.5)
    ax1.set_ylabel("dt [day]")
    ax1.set_title("Timestep Size vs Time")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    ax2.axvspan(shade_start, shade_end, color="0.85", alpha=0.5)
    ax2.plot(t_w, e_w, label="WHFast", lw=1.5)
    ax2.plot(t_i, e_i, label="IAS15", lw=1.5)
    ax2.plot(t_h, e_h, label="WHFast fake-HJ", lw=1.5)
    ax2.set_yscale("log")
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel("|relE|")
    ax2.set_title("Global Relative Energy Error vs Time")
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend()

    fig.tight_layout()
    fig.savefig(output_png, dpi=150)
    plt.close(fig)


def main() -> None:
    default_snapshot = SCRIPT_DIR / "find_close_encounter" / "apophis_1yr_before.bin"

    parser = argparse.ArgumentParser(
        description="Run WHFast and IAS15 from the saved Apophis snapshot and plot diagnostics."
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
        default=10.0,
        help="WHFast timestep in days.",
    )
    parser.add_argument("--whfast-csv", type=Path, default=SCRIPT_DIR / "whfast_python.csv")
    parser.add_argument("--ias15-csv", type=Path, default=SCRIPT_DIR / "ias15_python.csv")
    parser.add_argument("--fake-hj-csv", type=Path, default=SCRIPT_DIR / "whfast_fake_hj_python.csv")
    parser.add_argument("--plot", type=Path, default=SCRIPT_DIR / "integrator_compare_python.png")
    parser.add_argument(
        "--switch-start-yr",
        type=float,
        default=0.7,
        help="Time since restart to switch into fake hierarchical Jacobi ordering.",
    )
    parser.add_argument(
        "--switch-end-yr",
        type=float,
        default=1.1,
        help="Time since restart to switch back to the standard ordering.",
    )
    parser.add_argument(
        "--time-years",
        action="store_true",
        help="Plot the x-axis in years instead of days.",
    )
    args = parser.parse_args()

    snapshot_path = args.snapshot.resolve()
    if not snapshot_path.exists():
        raise FileNotFoundError(f"Snapshot file not found: {snapshot_path}")

    whfast_csv = args.whfast_csv.resolve()
    ias15_csv = args.ias15_csv.resolve()
    fake_hj_csv = args.fake_hj_csv.resolve()
    plot_path = args.plot.resolve()

    base_sim = rebound.Simulation(str(snapshot_path))
    rebound.clibrebound.reb_simulation_reset_function_pointers(byref(base_sim))
    for particle, key in zip(base_sim.particles[:3], STANDARD_ORDER):
        particle.hash = key

    print("Running WHFast from saved snapshot...", flush=True)
    whfast_summary = run_case(base_sim, "whfast", whfast_csv, args.duration_days, args.whfast_dt)
    print(
        f"Running WHFast fake-HJ window from {args.switch_start_yr:.3f} yr to {args.switch_end_yr:.3f} yr...",
        flush=True,
    )
    fake_hj_summary = run_fake_hj_window_case(
        base_sim,
        fake_hj_csv,
        args.duration_days,
        args.whfast_dt,
        args.switch_start_yr,
        args.switch_end_yr,
    )
    print("Running IAS15 from saved snapshot...", flush=True)
    ias15_summary = run_case(base_sim, "ias15", ias15_csv, args.duration_days, args.whfast_dt)
    print("Building comparison plot...", flush=True)
    make_plot(
        whfast_csv,
        ias15_csv,
        fake_hj_csv,
        plot_path,
        time_years=args.time_years,
        shade_start=switch_time_for_plot(fake_hj_summary, "switch_start_actual", args.time_years),
        shade_end=switch_time_for_plot(fake_hj_summary, "switch_end_actual", args.time_years),
    )

    print(f"Snapshot: {snapshot_path}")
    print(
        f"WHFast: steps={int(whfast_summary['steps_taken'])}, "
        f"t=[{whfast_summary['t_start']:.6f}, {whfast_summary['t_stop']:.6f}], "
        f"|relE_final|={whfast_summary['final_relE']:.6e}, "
        f"|relE|max={whfast_summary['max_relE']:.6e}"
    )
    print(
        f"WHFast fake-HJ: steps={int(fake_hj_summary['steps_taken'])}, "
        f"t=[{fake_hj_summary['t_start']:.6f}, {fake_hj_summary['t_stop']:.6f}], "
        f"|relE_final|={fake_hj_summary['final_relE']:.6e}, "
        f"|relE|max={fake_hj_summary['max_relE']:.6e}"
    )
    print(
        "WHFast fake-HJ actual window: "
        f"{(fake_hj_summary['switch_start_actual'] - fake_hj_summary['t_start']) / YEAR_DAYS:.6f} yr "
        f"to {(fake_hj_summary['switch_end_actual'] - fake_hj_summary['t_start']) / YEAR_DAYS:.6f} yr"
    )
    print(
        f"IAS15: steps={int(ias15_summary['steps_taken'])}, "
        f"t=[{ias15_summary['t_start']:.6f}, {ias15_summary['t_stop']:.6f}], "
        f"|relE_final|={ias15_summary['final_relE']:.6e}, "
        f"|relE|max={ias15_summary['max_relE']:.6e}"
    )
    print(f"Wrote {whfast_csv}")
    print(f"Wrote {fake_hj_csv}")
    print(f"Wrote {ias15_csv}")
    print(f"Wrote {plot_path}")


if __name__ == "__main__":
    main()
