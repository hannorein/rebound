#!/usr/bin/env python3
"""Resolved major-moon Solar System: WHFast versus fixed-tree WHFast-HJ.

This is a controlled Solar-System-like benchmark, not an ephemeris.  It uses
approximate present-day masses and orbital scales, deterministic orbital phases,
and resolves several dynamically important planet--moon branches:

    Earth--Moon
    Jupiter--Io--Europa--Ganymede--Callisto
    Saturn--Titan
    Uranus--Titania
    Neptune--Triton
    Pluto--Charon

The fixed HJ tree follows those physical branches.  Ordinary WHFast receives
the same Cartesian state and particle order, but its Jacobi coordinates form a
single comb and therefore cannot represent all resolved moon systems at once.
"""

from __future__ import annotations

import argparse
import csv
import inspect
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import time
import warnings


def find_repo_root() -> Path:
    for start in (Path.cwd().resolve(), Path(__file__).resolve().parent):
        for path in (start, *start.parents):
            if (path / "rebound").is_dir() and (path / "src").is_dir():
                return path
    raise RuntimeError("Could not locate the REBOUND repository root")


REPO_ROOT = find_repo_root()
sys.path.insert(0, str(REPO_ROOT))

CACHE = Path(tempfile.gettempdir()) / "rebound_complete_solar_cache"
CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(CACHE / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE))

import matplotlib  # noqa: E402
import numpy as np  # noqa: E402
import rebound  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


G = 4.0 * np.pi**2  # AU^3 / (solar mass yr^2)
KM_PER_AU = 149_597_870.7
DAYS_PER_YEAR = 365.25
ERROR_FLOOR = np.finfo(float).eps
OUTPUT_DIR = Path(__file__).resolve().parent


# Approximate masses in solar masses and heliocentric orbital elements.
# Angles are deliberately fixed rather than epoch-specific.
PLANETS = (
    ("Mercury", 1.6601e-7, 0.38710, 0.2056, 7.005, 15.0),
    ("Venus", 2.4478e-6, 0.72333, 0.0068, 3.395, 75.0),
    ("Earth", 3.0035e-6, 1.00000, 0.0167, 0.000, 130.0),
    ("Mars", 3.2272e-7, 1.52368, 0.0934, 1.850, 205.0),
    ("Jupiter", 9.5459e-4, 5.20260, 0.0489, 1.303, 250.0),
    ("Saturn", 2.8588e-4, 9.55490, 0.0565, 2.485, 310.0),
    ("Uranus", 4.3662e-5, 19.2184, 0.0463, 0.773, 25.0),
    ("Neptune", 5.1514e-5, 30.1104, 0.0095, 1.770, 95.0),
    ("Pluto", 6.55e-9, 39.482, 0.2488, 17.16, 170.0),
)


# name, mass [Msun], a [km], e, inclination to the planet's reference plane,
# deterministic true anomaly [deg].
MOONS = {
    "Earth": (("Moon", 3.694e-8, 384_400.0, 0.0549, 5.145, 35.0),),
    "Jupiter": (
        ("Io", 4.493e-8, 421_700.0, 0.0041, 0.04, 20.0),
        ("Europa", 2.413e-8, 671_100.0, 0.0094, 0.47, 110.0),
        ("Ganymede", 7.452e-8, 1_070_400.0, 0.0013, 0.20, 210.0),
        ("Callisto", 5.405e-8, 1_882_700.0, 0.0074, 0.19, 300.0),
    ),
    "Saturn": (("Titan", 6.764e-8, 1_221_870.0, 0.0288, 0.35, 70.0),),
    "Uranus": (("Titania", 1.77e-9, 435_910.0, 0.0011, 0.08, 145.0),),
    # Triton is retrograde; 156.9 degrees represents that geometry.
    "Neptune": (("Triton", 1.08e-8, 354_759.0, 0.00002, 156.9, 235.0),),
    "Pluto": (("Charon", 7.99e-10, 19_596.0, 0.0002, 0.0, 325.0),),
}


def detached_state(particle) -> np.ndarray:
    return np.asarray(
        [particle.x, particle.y, particle.z, particle.vx, particle.vy, particle.vz],
        dtype=float,
    )


def local_planet_system(planet_name: str, planet_mass: float):
    local = rebound.Simulation()
    local.G = G
    local.add(m=planet_mass, name=planet_name)
    for moon_name, mass, a_km, eccentricity, inclination, phase in MOONS.get(
        planet_name, ()
    ):
        local.add(
            m=mass,
            a=a_km / KM_PER_AU,
            e=eccentricity,
            inc=np.deg2rad(inclination),
            Omega=0.0,
            omega=0.0,
            f=np.deg2rad(phase),
            primary=local.particles[0],
            name=moon_name,
        )
    local.move_to_com()
    names = [planet_name] + [moon[0] for moon in MOONS.get(planet_name, ())]
    return names, [detached_state(p) for p in local.particles], sum(p.m for p in local.particles)


def nested_tree(labels: list[int]):
    tree: int | list = labels[0]
    for label in labels[1:]:
        tree = [tree, label]
    return tree


def compact_tree(tree) -> str:
    if isinstance(tree, int):
        return str(tree)
    return f"[{compact_tree(tree[0])},{compact_tree(tree[1])}]"


def initial_state():
    systems = {}
    for name, mass, *_ in PLANETS:
        systems[name] = local_planet_system(name, mass)

    # First place each planet-system barycentre on its heliocentric orbit.
    outer = rebound.Simulation()
    outer.G = G
    outer.add(m=1.0, name="Sun")
    for name, _, a, eccentricity, inclination, phase in PLANETS:
        outer.add(
            m=systems[name][2],
            a=a,
            e=eccentricity,
            inc=np.deg2rad(inclination),
            Omega=0.0,
            omega=0.0,
            f=np.deg2rad(phase),
            primary=outer.particles[0],
            name=f"{name}_system_com",
        )
    outer.move_to_com()

    names = ["Sun"]
    masses = [1.0]
    states = [detached_state(outer.particles[0])]
    branch_labels: dict[str, list[int]] = {}
    for outer_particle, (planet_name, planet_mass, *_rest) in zip(
        outer.particles[1:], PLANETS
    ):
        local_names, local_states, _ = systems[planet_name]
        labels = []
        center = detached_state(outer_particle)
        local_masses = [planet_mass] + [moon[1] for moon in MOONS.get(planet_name, ())]
        for local_name, mass, local_state in zip(local_names, local_masses, local_states):
            names.append(local_name)
            masses.append(mass)
            states.append(center + local_state)
            labels.append(len(names))  # HJ trees use one-based particle labels.
        branch_labels[planet_name] = labels

    # The solar hierarchy is a radial comb whose planet nodes may themselves be
    # branching planet--moon subtrees.
    solar_tree: int | list = 1
    for planet_name, *_ in PLANETS:
        labels = branch_labels[planet_name]
        branch = labels[0] if len(labels) == 1 else nested_tree(labels)
        solar_tree = [solar_tree, branch]

    return names, masses, states, compact_tree(solar_tree), branch_labels


NAMES, MASSES, STATES, HJ_TREE, BRANCH_LABELS = initial_state()


def build_simulation(mode: str, dt: float):
    sim = rebound.Simulation()
    sim.G = G
    for name, mass, values in zip(NAMES, MASSES, STATES):
        sim.add(
            m=mass,
            x=values[0],
            y=values[1],
            z=values[2],
            vx=values[3],
            vy=values[4],
            vz=values[5],
            name=name,
        )
    sim.move_to_com()
    sim.integrator = mode
    sim.dt = dt
    if mode == "whfast_hj":
        sim.integrate(sim.t, exact_finish_time=0, given_tree=True, tree=HJ_TREE)
    if mode == "ias15":
        sim.integrator.epsilon = 1.0e-12
    return sim


def state_array(sim) -> np.ndarray:
    return np.asarray([detached_state(p) for p in sim.particles])


def normalized_errors(candidate: np.ndarray, reference: np.ndarray):
    dx = np.sqrt(np.mean((candidate[:, :, :3] - reference[:, :, :3]) ** 2))
    dv = np.sqrt(np.mean((candidate[:, :, 3:] - reference[:, :, 3:]) ** 2))
    sx = max(np.sqrt(np.mean(reference[:, :, :3] ** 2)), ERROR_FLOOR)
    sv = max(np.sqrt(np.mean(reference[:, :, 3:] ** 2)), ERROR_FLOOR)
    return float(dx / sx), float(dv / sv)


def io_period() -> float:
    jupiter_mass = next(row[1] for row in PLANETS if row[0] == "Jupiter")
    io = MOONS["Jupiter"][0]
    return 2.0 * np.pi * np.sqrt((io[2] / KM_PER_AU) ** 3 / (G * (jupiter_mass + io[1])))


def run_fixed(mode: str, dt: float, duration: float, samples: int):
    sim = build_simulation(mode, dt)
    initial_energy = sim.energy()
    n_steps = int(math.ceil(duration / dt))
    stride = max(1, n_steps // max(1, samples - 1))
    times = [0.0]
    errors = [ERROR_FLOOR]
    states = [state_array(sim)]
    maximum = ERROR_FLOOR
    failure = None
    start = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            for step in range(1, n_steps + 1):
                sim.steps(1)
                error = abs((sim.energy() - initial_energy) / initial_energy)
                if not np.isfinite(error):
                    raise FloatingPointError("non-finite energy error")
                maximum = max(maximum, float(error))
                if step % stride == 0 or step == n_steps:
                    times.append(sim.t)
                    errors.append(max(float(error), ERROR_FLOOR))
                    states.append(state_array(sim))
        except Exception as exc:
            failure = f"{type(exc).__name__}: {exc}"
        if caught and failure is None:
            failure = f"{type(caught[0].message).__name__}: {caught[0].message}"
    return {
        "mode": mode,
        "dt": dt,
        "n_steps": n_steps,
        "times": np.asarray(times),
        "errors": np.asarray(errors),
        "states": np.asarray(states),
        "max_error": maximum,
        "runtime_seconds": time.perf_counter() - start,
        "failure": failure,
    }


def run_ias15(sample_times: np.ndarray):
    sim = build_simulation("ias15", sample_times[1] - sample_times[0])
    initial_energy = sim.energy()
    states = []
    errors = []
    start = time.perf_counter()
    for target in sample_times:
        sim.integrate(float(target), exact_finish_time=1)
        states.append(state_array(sim))
        errors.append(max(abs((sim.energy() - initial_energy) / initial_energy), ERROR_FLOOR))
    return {
        "mode": "ias15",
        "times": sample_times.copy(),
        "errors": np.asarray(errors),
        "states": np.asarray(states),
        "max_error": float(np.max(errors)),
        "runtime_seconds": time.perf_counter() - start,
        "failure": None,
    }


def save_main_plot(path: Path, whfast: dict, hj: dict, reference: dict):
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6))
    styles = (
        (whfast, "WHFast", "#d95f02"),
        (hj, "fixed-tree HJ", "#1b9e77"),
        (reference, "IAS15", "#3b5bdb"),
    )
    for run, label, color in styles:
        axes[0].semilogy(
            run["times"], run["errors"], color=color, lw=1.0, label=label
        )
    axes[0].set_xlabel("time [yr]")
    axes[0].set_ylabel(r"$|E(t)-E(0)|/|E(0)|$")
    axes[0].set_title("Total relative energy error")
    axes[0].legend()

    labels = ["WHFast", "fixed HJ"]
    values = [whfast["max_error"], hj["max_error"]]
    axes[1].bar(labels, values, color=["#d95f02", "#1b9e77"])
    axes[1].set_yscale("log")
    axes[1].set_ylabel("maximum relative energy error")
    axes[1].set_title(f"WHFast / HJ = {whfast['max_error']/hj['max_error']:.3g}x")
    for index, value in enumerate(values):
        axes[1].text(index, value * 1.15, f"{value:.3e}", ha="center", fontsize=9)
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
    fig.suptitle("Resolved major-moon Solar System")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_timestep_plot(path: Path, rows: list[dict]):
    x = np.asarray([1.0 / row["dt_divisor"] for row in rows])
    wh = np.asarray([row["whfast_max_error"] for row in rows])
    hj = np.asarray([row["hj_max_error"] for row in rows])
    order = np.argsort(x)
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4))
    axes[0].loglog(x[order], wh[order], "o-", color="#d95f02", label="WHFast")
    axes[0].loglog(x[order], hj[order], "s-", color="#1b9e77", label="fixed HJ")
    axes[0].set_xlabel(r"$\Delta t/P_{Io}$")
    axes[0].set_ylabel("maximum relative energy error")
    axes[0].legend()
    axes[1].loglog(x[order], (wh / hj)[order], "o-", color="#3b5bdb")
    axes[1].axhline(1.0, color="black", ls="--", lw=1)
    axes[1].set_xlabel(r"$\Delta t/P_{Io}$")
    axes[1].set_ylabel("WHFast error / fixed-HJ error")
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
    fig.suptitle("Timestep convergence: complete Solar benchmark")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--years", type=float, default=5.0)
    parser.add_argument("--divisors", type=int, nargs="+", default=(10, 20, 40, 80))
    parser.add_argument("--main-divisor", type=int, default=20)
    parser.add_argument("--samples", type=int, default=1201)
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.years <= 0 or args.samples < 2 or any(value <= 0 for value in args.divisors):
        raise ValueError("years/divisors must be positive and samples at least two")
    if args.main_divisor not in args.divisors:
        raise ValueError("main divisor must be present in divisors")
    if "given_tree" not in inspect.signature(rebound.Simulation.integrate).parameters:
        raise RuntimeError("The loaded REBOUND build lacks the fixed-tree HJ API")

    args.output.mkdir(parents=True, exist_ok=True)
    period = io_period()
    rows = []
    main_runs = None
    print(f"particles={len(NAMES)}, P_Io={period:.9g} yr ({period*DAYS_PER_YEAR:.6g} d)")
    print("tree:", HJ_TREE)
    for divisor in args.divisors:
        dt = period / divisor
        whfast = run_fixed("whfast", dt, args.years, args.samples)
        hj = run_fixed("whfast_hj", dt, args.years, args.samples)
        if whfast["failure"] or hj["failure"]:
            raise RuntimeError(f"P/{divisor} failed: {whfast['failure'] or hj['failure']}")
        row = {
            "dt_divisor": divisor,
            "dt_years": dt,
            "dt_days": dt * DAYS_PER_YEAR,
            "years": args.years,
            "n_steps": whfast["n_steps"],
            "whfast_max_error": whfast["max_error"],
            "hj_max_error": hj["max_error"],
            "improvement": whfast["max_error"] / hj["max_error"],
            "whfast_runtime_seconds": whfast["runtime_seconds"],
            "hj_runtime_seconds": hj["runtime_seconds"],
        }
        rows.append(row)
        print(
            f"P_Io/{divisor:<3d}: WH={row['whfast_max_error']:.3e} "
            f"HJ={row['hj_max_error']:.3e} ratio={row['improvement']:.3e}"
        )
        if divisor == args.main_divisor:
            main_runs = (whfast, hj)

    assert main_runs is not None
    whfast, hj = main_runs
    reference = run_ias15(whfast["times"])
    wh_dx, wh_dv = normalized_errors(whfast["states"], reference["states"])
    hj_dx, hj_dv = normalized_errors(hj["states"], reference["states"])
    selected = dict(next(row for row in rows if row["dt_divisor"] == args.main_divisor))
    selected.update(
        {
            "ias15_max_error": reference["max_error"],
            "whfast_position_error": wh_dx,
            "whfast_velocity_error": wh_dv,
            "hj_position_error": hj_dx,
            "hj_velocity_error": hj_dv,
            "ias15_runtime_seconds": reference["runtime_seconds"],
        }
    )

    with (args.output / "timestep_results.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (args.output / "main_result.json").open("w") as stream:
        json.dump(selected, stream, indent=2)
    with (args.output / "setup.json").open("w") as stream:
        json.dump(
            {
                "G": G,
                "units": {"length": "AU", "mass": "solar mass", "time": "year"},
                "names": NAMES,
                "masses": MASSES,
                "tree": HJ_TREE,
                "planet_data": PLANETS,
                "moon_data": MOONS,
                "io_period_years": period,
            },
            stream,
            indent=2,
        )

    with (args.output / "energy_timeseries.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("method", "time_years", "relative_energy_error"))
        for run in (whfast, hj, reference):
            for t, error in zip(run["times"], run["errors"]):
                writer.writerow((run["mode"], t, error))

    save_main_plot(args.output / "energy_comparison.png", whfast, hj, reference)
    save_timestep_plot(args.output / "timestep_convergence.png", rows)
    print("IAS15 max error:", f"{reference['max_error']:.3e}")
    print(f"trajectory dx: WH={wh_dx:.3e}, HJ={hj_dx:.3e}")
    print(f"trajectory dv: WH={wh_dv:.3e}, HJ={hj_dv:.3e}")
    print("outputs:", args.output)


if __name__ == "__main__":
    main()
