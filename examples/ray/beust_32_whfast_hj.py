#!/usr/bin/env python3
"""Reproduce Beust section 3.2 triple-star tests with whfast_hj.

This script keeps the Hierarchical Jacobi tree fixed as [[1,2],3].
It reports hierarchy diagnostics, but never uses them to choose a new tree.
"""

from pathlib import Path
import argparse
import csv
import inspect
import os
import sys
import tempfile


def find_repo_root():
    starts = [Path.cwd().resolve(), Path(__file__).resolve().parent]
    seen = set()
    for start in starts:
        for path in (start, *start.parents):
            if path in seen:
                continue
            seen.add(path)
            if (path / "rebound").is_dir() and (path / "src").is_dir():
                return path
    raise RuntimeError("Could not find the REBOUND repository root.")


repo_root = find_repo_root()
sys.path.insert(0, str(repo_root))

import rebound  # noqa: E402
import ctypes  # noqa: E402,F401

_cache_root = Path(tempfile.gettempdir()) / "rebound_beust_32_cache"
_cache_root.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", str(_cache_root))
os.environ.setdefault(
    "MPLCONFIGDIR",
    str(_cache_root / "matplotlib"),
)

import matplotlib  # noqa: E402
import numpy as np  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


HJ_TREE = "[[1,2],3]"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "beust_32_outputs"

G = 1.0
A_INNER = 1.0
A_OUTER = 10.0
E_INNER = 0.08
E_OUTER = 0.27
Q_INNER = 0.27
Q_OUTER = 0.44
OUTER_PERIASTRON_OFFSET_DEG = 270.0


def require_given_tree_api():
    integrate_signature = inspect.signature(rebound.Simulation.integrate)
    if "given_tree" not in integrate_signature.parameters:
        raise RuntimeError(
            "This process loaded an old rebound module without given_tree. "
            "Check PYTHONPATH or rebuild/reinstall the local REBOUND package."
        )
    print("rebound python:", rebound.__file__)
    print("rebound library:", rebound.__libpath__)
    print("integrate:", integrate_signature)


def inner_period():
    return 2.0 * np.pi * np.sqrt(A_INNER**3 / (G * 1.0))


def actual_outer_period_ratio():
    p_inner = inner_period()
    p_outer = 2.0 * np.pi * np.sqrt(A_OUTER**3 / (G * (1.0 + Q_OUTER)))
    return p_outer / p_inner


def setup_beust_triple(mutual_inclination_deg):
    """Set up Beust Sect. 3.2's synthetic hierarchical triple."""
    m2 = 1.0 / (1.0 + Q_INNER)
    m1 = Q_INNER * m2
    m3 = Q_OUTER * (m1 + m2)

    sim = rebound.Simulation()
    sim.G = G
    sim.add(m=m1, name="inner_1")
    sim.add(
        m=m2,
        a=A_INNER,
        e=E_INNER,
        inc=0.0,
        Omega=0.0,
        omega=0.0,
        f=0.0,
        primary=sim.particles[0],
        name="inner_2",
    )

    inner_com = sim.com(last=2)
    sim.add(
        m=m3,
        a=A_OUTER,
        e=E_OUTER,
        inc=np.deg2rad(mutual_inclination_deg),
        Omega=0.0,
        omega=np.deg2rad(OUTER_PERIASTRON_OFFSET_DEG),
        f=0.0,
        primary=inner_com,
        name="outer_3",
    )

    sim.move_to_com()
    sim.integrator = "whfast_hj"
    sim.dt = inner_period() / 20.0
    return sim


def vector_from_particle(particle, origin, attr):
    return np.array(
        [
            getattr(particle, attr[0]) - getattr(origin, attr[0]),
            getattr(particle, attr[1]) - getattr(origin, attr[1]),
            getattr(particle, attr[2]) - getattr(origin, attr[2]),
        ],
        dtype=float,
    )


def relative_state(particle, origin):
    r = vector_from_particle(particle, origin, ("x", "y", "z"))
    v = vector_from_particle(particle, origin, ("vx", "vy", "vz"))
    return r, v


def mutual_inclination_deg(sim):
    particles = sim.particles
    inner_r, inner_v = relative_state(particles[1], particles[0])
    outer_r, outer_v = relative_state(particles[2], sim.com(last=2))

    h_inner = np.cross(inner_r, inner_v)
    h_outer = np.cross(outer_r, outer_v)
    denom = np.linalg.norm(h_inner) * np.linalg.norm(h_outer)
    if denom == 0.0:
        return np.nan

    cos_i = np.dot(h_inner, h_outer) / denom
    cos_i = np.clip(cos_i, -1.0, 1.0)
    return float(np.rad2deg(np.arccos(cos_i)))


def orbital_diagnostics(sim):
    particles = sim.particles
    inner_orbit = particles[1].orbit(primary=particles[0])
    outer_orbit = particles[2].orbit(primary=sim.com(last=2))

    inner_apocenter = inner_orbit.a * (1.0 + inner_orbit.e)
    outer_pericenter = outer_orbit.a * (1.0 - outer_orbit.e)
    if inner_apocenter > 0.0:
        hierarchy_ratio = outer_pericenter / inner_apocenter
    else:
        hierarchy_ratio = np.nan

    return {
        "e_inner": inner_orbit.e,
        "mutual_inclination_deg": mutual_inclination_deg(sim),
        "a_inner": inner_orbit.a,
        "a_outer": outer_orbit.a,
        "hierarchy_ratio": hierarchy_ratio,
    }


def run_case(mutual_inclination, n_periods, n_outputs):
    sim = setup_beust_triple(mutual_inclination)
    period = inner_period()
    energy_initial = sim.energy()
    output_times = np.linspace(0.0, n_periods * period, n_outputs)

    rows = []
    for target_time in output_times:
        if target_time > sim.t:
            sim.integrate(target_time, exact_finish_time=0, given_tree=True, tree=HJ_TREE)

        diagnostics = orbital_diagnostics(sim)
        energy = sim.energy()
        rel_energy_error = abs((energy - energy_initial) / energy_initial)
        rows.append(
            {
                "time_inner_periods": sim.t / period,
                "e_inner": diagnostics["e_inner"],
                "mutual_inclination_deg": diagnostics["mutual_inclination_deg"],
                "rel_energy_error": rel_energy_error,
                "a_inner": diagnostics["a_inner"],
                "a_outer": diagnostics["a_outer"],
                "hierarchy_ratio": diagnostics["hierarchy_ratio"],
            }
        )

    return rows


def write_csv(path, rows):
    fieldnames = [
        "time_inner_periods",
        "e_inner",
        "mutual_inclination_deg",
        "rel_energy_error",
        "a_inner",
        "a_outer",
        "hierarchy_ratio",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def column(rows, name):
    return np.array([row[name] for row in rows], dtype=float)


def save_plot(path, rows, y_name, ylabel, title, semilogy=False):
    t = column(rows, "time_inner_periods")
    y = column(rows, y_name)
    if semilogy:
        y = np.where(y > 0.0, y, np.nan)

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    if semilogy:
        ax.semilogy(t, y, color="black", linewidth=0.8)
    else:
        ax.plot(t, y, color="black", linewidth=0.8)
    ax.set_xlabel("time (inner revolutions)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def summarize_case(mutual_inclination, rows):
    e_inner = column(rows, "e_inner")
    mutual_i = column(rows, "mutual_inclination_deg")
    rel_e = column(rows, "rel_energy_error")
    hierarchy = column(rows, "hierarchy_ratio")
    print(
        f"i={mutual_inclination:g} deg: "
        f"e_inner=[{np.nanmin(e_inner):.4g}, {np.nanmax(e_inner):.4g}], "
        f"mutual_i=[{np.nanmin(mutual_i):.3g}, {np.nanmax(mutual_i):.3g}] deg, "
        f"max dE/E={np.nanmax(rel_e):.3e}, "
        f"min hierarchy ratio={np.nanmin(hierarchy):.3g}"
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run Beust section 3.2 fixed-tree whfast_hj triple tests."
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--quick",
        action="store_true",
        help="Run a short validation integration. This is the default.",
    )
    mode.add_argument(
        "--paper",
        action="store_true",
        help="Run the Beust-style 20000-inner-period integration.",
    )
    parser.add_argument(
        "--n-periods",
        type=float,
        default=None,
        help="Override the integration length in inner binary periods.",
    )
    parser.add_argument(
        "--n-outputs",
        type=int,
        default=None,
        help="Override the number of saved output samples.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for CSV and PNG outputs.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    require_given_tree_api()

    if args.n_periods is None:
        n_periods = 20000.0 if args.paper else 500.0
    else:
        n_periods = args.n_periods

    if args.n_outputs is None:
        n_outputs = 2001 if args.paper else 501
    else:
        n_outputs = args.n_outputs
    if n_outputs < 2:
        raise ValueError("--n-outputs must be at least 2.")

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print("fixed whfast_hj tree:", HJ_TREE)
    print(f"dt = P_inner/20 = {inner_period()/20.0:.16g}")
    print(f"a_inner/a_outer = {A_INNER/A_OUTER:.6g}")
    print(f"actual P_outer/P_inner = {actual_outer_period_ratio():.6g}")
    print(f"n_periods = {n_periods:g}, n_outputs = {n_outputs}")
    print("output dir:", output_dir)

    results = {}
    for mutual_inclination in (20.0, 60.0):
        rows = run_case(mutual_inclination, n_periods, n_outputs)
        results[mutual_inclination] = rows
        csv_path = output_dir / f"beust_32_i{int(mutual_inclination):02d}.csv"
        write_csv(csv_path, rows)
        summarize_case(mutual_inclination, rows)
        print("wrote", csv_path)

    save_plot(
        output_dir / "beust_32_i20_e_inner.png",
        results[20.0],
        "e_inner",
        "e_inner",
        "Beust 3.2 fixed-tree whfast_hj, i = 20 deg",
    )
    save_plot(
        output_dir / "beust_32_i20_energy_error.png",
        results[20.0],
        "rel_energy_error",
        "abs((E - E0) / E0)",
        "Beust 3.2 fixed-tree whfast_hj energy error, i = 20 deg",
        semilogy=True,
    )
    save_plot(
        output_dir / "beust_32_i60_e_inner.png",
        results[60.0],
        "e_inner",
        "e_inner",
        "Beust 3.2 fixed-tree whfast_hj, i = 60 deg",
    )
    save_plot(
        output_dir / "beust_32_i60_mutual_inclination.png",
        results[60.0],
        "mutual_inclination_deg",
        "mutual inclination (deg)",
        "Beust 3.2 fixed-tree whfast_hj, i = 60 deg",
    )
    print("wrote plots to", output_dir)


if __name__ == "__main__":
    main()
