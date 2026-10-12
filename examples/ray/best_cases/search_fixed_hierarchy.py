#!/usr/bin/env python3
"""Search for stable fixed hierarchies where WHFast-HJ beats WHFast.

The primary family is a 2+2 stellar quadruple.  A standard Jacobi comb can
make the first tight binary a Kepler problem, but not the second one.  The
fixed HJ tree ``[[1,2],[3,4]]`` makes both binaries (and their wide mutual
orbit) Kepler problems.  Both methods use the same physical initial state,
fixed timestep, second-order map, integration duration, and output cadence.

Outputs are kept under ``best_5/``.  Each ranked case contains its parameters,
the two energy-error time series, and a plot.  ``ranking.csv`` contains all
successful candidates, so a longer run can replace the current top five.
"""

from __future__ import annotations

import argparse
import csv
import inspect
import itertools
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import time
import warnings


def find_repo_root() -> Path:
    starts = (Path.cwd().resolve(), Path(__file__).resolve().parent)
    for start in starts:
        for path in (start, *start.parents):
            if (path / "rebound").is_dir() and (path / "src").is_dir():
                return path
    raise RuntimeError("Could not locate the REBOUND repository root")


REPO_ROOT = find_repo_root()
sys.path.insert(0, str(REPO_ROOT))

cache = Path(tempfile.gettempdir()) / "rebound_hj_energy_search_cache"
cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(cache))

import matplotlib  # noqa: E402
import numpy as np  # noqa: E402
import rebound  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


TREE = "[[1,2],[3,4]]"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "best_5"
ERROR_FLOOR = np.finfo(float).eps


def kepler_period(a: float, mass: float) -> float:
    return 2.0 * np.pi * np.sqrt(a**3 / mass)


def binary_offsets(m1, m2, a, e, inc, omega, f):
    local = rebound.Simulation()
    local.G = 1.0
    local.add(m=m1)
    local.add(
        m=m2,
        a=a,
        e=e,
        inc=inc,
        Omega=0.0,
        omega=omega,
        f=f,
        primary=local.particles[0],
    )
    local.move_to_com()
    return [
        np.array([p.x, p.y, p.z, p.vx, p.vy, p.vz], dtype=float)
        for p in local.particles
    ]


def build_simulation(config: dict, mode: str) -> rebound.Simulation:
    masses = config["masses"]
    a1, a2 = config["a_inner"]
    e1, e2 = config["e_inner"]
    phase1, phase2, outer_phase = np.deg2rad(config["phases_deg"])
    mutual_inc = np.deg2rad(config["mutual_inc_deg"])

    left = binary_offsets(masses[0], masses[1], a1, e1, 0.0, 0.0, phase1)
    right = binary_offsets(
        masses[2], masses[3], a2, e2, mutual_inc, 0.37, phase2
    )
    outer = binary_offsets(
        masses[0] + masses[1],
        masses[2] + masses[3],
        config["a_outer"],
        config["e_outer"],
        0.5 * mutual_inc,
        0.71,
        outer_phase,
    )

    sim = rebound.Simulation()
    sim.G = 1.0
    for state, center, mass, name in zip(
        left + right,
        [outer[0], outer[0], outer[1], outer[1]],
        masses,
        ("A1", "A2", "B1", "B2"),
    ):
        combined = state + center
        sim.add(
            m=mass,
            x=combined[0],
            y=combined[1],
            z=combined[2],
            vx=combined[3],
            vy=combined[4],
            vz=combined[5],
            name=name,
        )
    sim.move_to_com()
    sim.integrator = mode
    sim.dt = config["dt"]
    if mode == "whfast_hj":
        # Parse and retain the supplied fixed tree before the first real step.
        sim.integrate(sim.t, exact_finish_time=0, given_tree=True, tree=TREE)
    return sim


def hierarchy_margin(config: dict) -> float:
    inner_apocenter = max(
        a * (1.0 + e) for a, e in zip(config["a_inner"], config["e_inner"])
    )
    return config["a_outer"] * (1.0 - config["e_outer"]) / inner_apocenter


def run_method(config: dict, mode: str, samples: int):
    sim = build_simulation(config, mode)
    initial_energy = sim.energy()
    n_steps = config["n_steps"]
    sample_stride = max(1, n_steps // max(1, samples - 1))
    times = [0.0]
    errors = [ERROR_FLOOR]
    maximum_error = ERROR_FLOOR
    failed = None
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            for step in range(1, n_steps + 1):
                sim.steps(1)
                error = abs((sim.energy() - initial_energy) / initial_energy)
                if not np.isfinite(error):
                    raise FloatingPointError("non-finite relative energy error")
                error = max(float(error), ERROR_FLOOR)
                maximum_error = max(maximum_error, error)
                if step % sample_stride == 0 or step == n_steps:
                    times.append(sim.t)
                    errors.append(error)
        except Exception as exc:  # retain failed parameter sets in the audit
            failed = f"{type(exc).__name__}: {exc}"
        if caught and failed is None:
            failed = f"{type(caught[0].message).__name__}: {caught[0].message}"
    return np.asarray(times), np.asarray(errors), maximum_error, failed


def make_config(index, masses, a2, e_pair, outer_a, outer_e, inc, divisor, phases):
    a1 = 1.0
    p1 = kepler_period(a1, masses[0] + masses[1])
    p2 = kepler_period(a2, masses[2] + masses[3])
    reference_period = min(p1, p2)
    duration = 120.0 * max(p1, p2)
    dt = reference_period / divisor
    return {
        "id": f"candidate_{index:04d}",
        "family": "stable_2plus2",
        "tree": TREE,
        "masses": list(masses),
        "a_inner": [a1, a2],
        "e_inner": list(e_pair),
        "a_outer": outer_a,
        "e_outer": outer_e,
        "mutual_inc_deg": inc,
        "phases_deg": list(phases),
        "dt_divisor": divisor,
        "reference_period": reference_period,
        "dt": dt,
        "duration": duration,
        "n_steps": int(math.ceil(duration / dt)),
    }


def candidate_configs(limit: int):
    mass_sets = (
        (1.0, 1.0, 1.0, 1.0),
        (1.0, 0.7, 0.9, 0.6),
        (1.0, 0.2, 0.8, 0.15),
        (1.3, 0.3, 0.7, 0.7),
    )
    a2_values = (0.55, 0.8, 1.0, 1.35)
    eccentricities = ((0.0, 0.0), (0.1, 0.2), (0.35, 0.1), (0.5, 0.4))
    outer_orbits = ((12.0, 0.0), (20.0, 0.2), (40.0, 0.35), (80.0, 0.5))
    inclinations = (0.0, 35.0, 80.0, 145.0)
    divisors = (10, 20, 40)
    phases = ((0.0, 137.0, 61.0), (73.0, 211.0, 179.0))

    # Interleave the grid deterministically so even short searches cover every
    # physical axis rather than exhausting the first mass set.
    raw = list(
        itertools.product(
            mass_sets,
            a2_values,
            eccentricities,
            outer_orbits,
            inclinations,
            divisors,
            phases,
        )
    )
    rng = np.random.default_rng(20260810)
    rng.shuffle(raw)
    for index, values in enumerate(raw[:limit], 1):
        masses, a2, e_pair, outer_orbit, inc, divisor, phase_tuple = values
        config = make_config(
            index,
            masses,
            a2,
            e_pair,
            outer_orbit[0],
            outer_orbit[1],
            inc,
            divisor,
            phase_tuple,
        )
        if hierarchy_margin(config) >= 6.0:
            yield config


def evaluate(config: dict, samples: int) -> dict:
    start = time.perf_counter()
    result = {**config, "hierarchy_margin": hierarchy_margin(config)}
    series = {}
    for mode in ("whfast", "whfast_hj"):
        times, errors, maximum_error, failed = run_method(config, mode, samples)
        series[mode] = (times, errors)
        result[f"{mode}_failure"] = failed
        result[f"{mode}_max_error"] = maximum_error if failed is None else np.nan
        result[f"{mode}_rms_error"] = (
            float(np.sqrt(np.mean(errors**2))) if failed is None else np.nan
        )
    if any(result[f"{mode}_failure"] for mode in ("whfast", "whfast_hj")):
        result["improvement"] = np.nan
    else:
        result["improvement"] = (
            result["whfast_max_error"] / result["whfast_hj_max_error"]
        )
    result["runtime_seconds"] = time.perf_counter() - start
    result["series"] = series
    return result


def scalar_result(result: dict) -> dict:
    return {k: v for k, v in result.items() if k != "series"}


def write_ranking(path: Path, results: list[dict]):
    rows = [scalar_result(r) for r in results]
    fieldnames = list(rows[0]) if rows else []
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_best_case(rank: int, result: dict, output: Path):
    case_dir = output / f"rank_{rank:02d}"
    case_dir.mkdir(parents=True, exist_ok=True)
    # A new search can change which physical case occupies a rank. Remove only
    # derived audits whose metadata would otherwise refer to the previous case.
    for stale_name in (
        "timestep_validation.png",
        "whfast_order_audit.csv",
        "ias15_validation.json",
        "ias15_validation.png",
    ):
        stale_path = case_dir / stale_name
        if stale_path.exists():
            stale_path.unlink()
    with (case_dir / "parameters.json").open("w") as stream:
        json.dump(scalar_result(result), stream, indent=2)

    modes = ("whfast", "whfast_hj")
    series_by_mode = result["series"]
    with (case_dir / "energy_errors.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("method", "time", "time_over_reference_period", "rel_energy_error"))
        for mode in modes:
            times, errors = series_by_mode[mode]
            for t, error in zip(times, errors):
                writer.writerow((mode, t, t / result["reference_period"], error))

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    labels = {"whfast": "WHFast (Jacobi comb)", "whfast_hj": "WHFast-HJ fixed tree"}
    colors = {"whfast": "#d95f02", "whfast_hj": "#1b9e77"}
    for mode in modes:
        times, errors = series_by_mode[mode]
        ax.semilogy(
            times / result["reference_period"],
            errors,
            lw=1.0,
            color=colors[mode],
            label=labels[mode],
        )
    ax.set_xlabel(r"time / $P_{\rm ref}$")
    ax.set_ylabel(r"$|E(t)-E(0)|/|E(0)|$")
    ax.set_title(
        f"Rank {rank}: fixed HJ improves max energy error by "
        f"{result['improvement']:.3g}x"
    )
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(case_dir / "energy_error.png", dpi=180)
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=int, default=160)
    parser.add_argument("--samples", type=int, default=1201)
    parser.add_argument("--top", type=int, default=5)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--stop-after-target",
        action="store_true",
        help="Stop after finding five cases above the 10,000x target.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.candidates < 1 or args.samples < 2 or args.top < 1:
        raise ValueError("candidates and top must be positive; samples must be at least 2")
    if "given_tree" not in inspect.signature(rebound.Simulation.integrate).parameters:
        raise RuntimeError("The loaded REBOUND build does not provide the fixed-tree HJ API")

    args.output.mkdir(parents=True, exist_ok=True)
    results = []
    target_count = 0
    print("REBOUND:", rebound.__file__)
    print("fixed HJ tree:", TREE)
    for number, config in enumerate(candidate_configs(args.candidates), 1):
        result = evaluate(config, args.samples)
        results.append(result)
        improvement = result["improvement"]
        if np.isfinite(improvement) and improvement >= 1.0e4:
            target_count += 1
        print(
            f"[{number:03d}] {config['id']} margin={result['hierarchy_margin']:.2f} "
            f"dt=P/{config['dt_divisor']:d} WH={result['whfast_max_error']:.3e} "
            f"HJ={result['whfast_hj_max_error']:.3e} ratio={improvement:.3e}"
        )
        if args.stop_after_target and target_count >= args.top:
            break

    results.sort(
        key=lambda row: row["improvement"] if np.isfinite(row["improvement"]) else -np.inf,
        reverse=True,
    )
    write_ranking(args.output / "ranking.csv", results)
    best = [r for r in results if np.isfinite(r["improvement"])][: args.top]
    for rank, result in enumerate(best, 1):
        write_best_case(rank, result, args.output)

    summary = {
        "target_improvement": 1.0e4,
        "evaluated": len(results),
        "target_count": sum(r["improvement"] >= 1.0e4 for r in best),
        "top_cases": [scalar_result(r) for r in best],
    }
    with (args.output / "summary.json").open("w") as stream:
        json.dump(summary, stream, indent=2)
    print("\nTop cases:")
    for rank, result in enumerate(best, 1):
        print(
            f"  {rank}. {result['id']}: {result['improvement']:.6g}x "
            f"(WH={result['whfast_max_error']:.3e}, "
            f"HJ={result['whfast_hj_max_error']:.3e})"
        )
    print("outputs:", args.output)


if __name__ == "__main__":
    main()
