#!/usr/bin/env python3
"""Qualitative Beust section 3.3 GG Tau test for whfast_hj.

The experiment follows a fixed four-star hierarchy,

    [[1,2],[3,4]]

and attaches massless circumbinary disk particles to the central binary
side of the tree. The tree is fixed during the integration; diagnostics
never trigger hierarchy changes.
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

_cache_root = Path(tempfile.gettempdir()) / "rebound_beust_33_cache"
_cache_root.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", str(_cache_root))
os.environ.setdefault("MPLCONFIGDIR", str(_cache_root / "matplotlib"))

import rebound  # noqa: E402
import matplotlib  # noqa: E402
import numpy as np  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "beust_33_outputs"
G_AU_YR_MSUN = 4.0 * np.pi * np.pi

M_CENTRAL_1 = 0.6
M_CENTRAL_2 = 0.7
M_OUTER_1 = 0.12
M_OUTER_2 = 0.04

CENTRAL_A_AU = 35.0
CENTRAL_E = 0.30
OUTER_BINARY_A_AU = 190.0
OUTER_BINARY_E = 0.20

WIDE_PERIASTRON_AU = 833.0
WIDE_E = 0.50
WIDE_I_DEG = 20.0
OUTER_BINARY_REL_I_DEG = 150.0
FIRST_WIDE_PERIASTRON_YR = 44000.0

DISK_R_MIN_AU = 150.0
DISK_R_MAX_AU = 1200.0
DISK_E_MAX = 0.05
DISK_INC_MAX_DEG = 2.0

MAX_EXPLICIT_DISK_PARTICLES = 5000


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


def kepler_period(a, total_mass):
    return np.sqrt(a**3 / total_mass)


def wide_semimajor_axis(wide_periastron, wide_e):
    if not (0.0 <= wide_e < 1.0):
        raise ValueError("wide eccentricity must satisfy 0 <= e < 1.")
    return wide_periastron / (1.0 - wide_e)


def mean_anomaly_for_next_periastron(period, time_to_periastron):
    phase_to_periastron = (time_to_periastron / period) % 1.0
    return 2.0 * np.pi * (1.0 - phase_to_periastron)


def particle_state(particle):
    return np.array([particle.x, particle.y, particle.z]), np.array(
        [particle.vx, particle.vy, particle.vz]
    )


def add_outer_binary_at_wide_com(
    sim,
    wide_com,
    outer_binary_a,
    outer_binary_e,
    wide_i_deg,
    outer_binary_rel_i_deg,
):
    local = rebound.Simulation()
    local.G = sim.G
    local.add(m=M_OUTER_1, name="outer_1_local")
    local.add(
        m=M_OUTER_2,
        a=outer_binary_a,
        e=outer_binary_e,
        inc=np.deg2rad(wide_i_deg + outer_binary_rel_i_deg),
        Omega=0.0,
        omega=0.0,
        f=0.0,
        primary=local.particles[0],
        name="outer_2_local",
    )
    local.move_to_com()

    for p_local, name in zip(local.particles, ("outer_1", "outer_2")):
        sim.add(
            m=p_local.m,
            x=wide_com.x + p_local.x,
            y=wide_com.y + p_local.y,
            z=wide_com.z + p_local.z,
            vx=wide_com.vx + p_local.vx,
            vy=wide_com.vy + p_local.vy,
            vz=wide_com.vz + p_local.vz,
            name=name,
        )


def sample_disk_elements(rng, n_disk):
    r = rng.uniform(DISK_R_MIN_AU, DISK_R_MAX_AU, n_disk)
    e = rng.uniform(0.0, DISK_E_MAX, n_disk)
    inc = rng.uniform(0.0, np.deg2rad(DISK_INC_MAX_DEG), n_disk)
    Omega = rng.uniform(0.0, 2.0 * np.pi, n_disk)
    omega = rng.uniform(0.0, 2.0 * np.pi, n_disk)
    f = rng.uniform(0.0, 2.0 * np.pi, n_disk)
    a = r * (1.0 + e * np.cos(f)) / (1.0 - e * e)
    return a, e, inc, Omega, omega, f


def add_circumbinary_disk(sim, n_disk, seed):
    rng = np.random.default_rng(seed)
    central_com = sim.com(last=2)
    elements = sample_disk_elements(rng, n_disk)
    for a, e, inc, Omega, omega, f in zip(*elements):
        sim.add(
            m=0.0,
            a=float(a),
            e=float(e),
            inc=float(inc),
            Omega=float(Omega),
            omega=float(omega),
            f=float(f),
            primary=central_com,
        )


def setup_ggtau_sim(args):
    sim = rebound.Simulation()
    sim.G = G_AU_YR_MSUN

    sim.add(m=M_CENTRAL_1, name="central_1")
    sim.add(
        m=M_CENTRAL_2,
        a=args.central_a,
        e=args.central_e,
        inc=0.0,
        Omega=0.0,
        omega=0.0,
        f=0.0,
        primary=sim.particles[0],
        name="central_2",
    )

    central_com = sim.com(last=2)
    wide_a = wide_semimajor_axis(args.wide_periastron, args.wide_e)
    total_outer_mass = M_OUTER_1 + M_OUTER_2
    wide_period = kepler_period(wide_a, M_CENTRAL_1 + M_CENTRAL_2 + total_outer_mass)
    wide_M = mean_anomaly_for_next_periastron(
        wide_period, args.first_periastron_time
    )
    wide_com = rebound.Particle(
        simulation=sim,
        primary=central_com,
        m=total_outer_mass,
        a=wide_a,
        e=args.wide_e,
        inc=np.deg2rad(args.wide_i),
        Omega=0.0,
        omega=0.0,
        M=wide_M,
    )
    add_outer_binary_at_wide_com(
        sim,
        wide_com,
        args.outer_binary_a,
        args.outer_binary_e,
        args.wide_i,
        args.outer_binary_rel_i,
    )

    add_circumbinary_disk(sim, args.n_disk, args.seed)
    sim.N_active = 4
    sim.testparticle_type = 0
    sim.move_to_com()
    sim.integrator = "whfast_hj"
    sim.dt = kepler_period(args.central_a, M_CENTRAL_1 + M_CENTRAL_2) / args.dt_divisor
    return sim


def make_fixed_tree(n_particles):
    central_side = "[1,2]"
    for label in range(5, n_particles + 1):
        central_side = f"[{central_side},{label}]"
    return f"[{central_side},[3,4]]"


def relative_positions(sim, first, last, origin):
    particles = sim.particles
    rows = []
    for i in range(first, last):
        p = particles[i]
        rows.append((p.x - origin.x, p.y - origin.y, p.z - origin.z))
    if not rows:
        return np.zeros((0, 3))
    return np.array(rows, dtype=float)


def orbit_or_nan(particle, primary, G=None):
    try:
        if G is None:
            return particle.orbit(primary=primary)
        return particle.orbit(primary=primary, G=G)
    except Exception:
        return None


def disk_radii(sim):
    central_com = sim.com(last=2)
    xyz = relative_positions(sim, 4, sim.N, central_com)
    return np.sqrt(np.sum(xyz * xyz, axis=1)), xyz


def disk_edge_summary(radii, survivor_radius):
    finite = radii[np.isfinite(radii)]
    within = finite[finite <= survivor_radius]
    if within.size == 0:
        return {
            "n_within_survivor_radius": 0,
            "inner_edge_p05": np.nan,
            "outer_edge_p95": np.nan,
            "median_radius": np.nan,
        }
    return {
        "n_within_survivor_radius": int(within.size),
        "inner_edge_p05": float(np.percentile(within, 5)),
        "outer_edge_p95": float(np.percentile(within, 95)),
        "median_radius": float(np.median(within)),
    }


def snapshot_summary(sim, energy_initial, survivor_radius):
    ps = sim.particles
    central_com = sim.com(last=2)
    outer_com = sim.com(first=2, last=4)
    central_orbit = orbit_or_nan(ps[1], ps[0])
    outer_binary_orbit = orbit_or_nan(ps[3], ps[2])
    wide_orbit = orbit_or_nan(outer_com, central_com, G=sim.G)
    radii, _ = disk_radii(sim)
    disk_summary = disk_edge_summary(radii, survivor_radius)

    def value(orbit, attr):
        return np.nan if orbit is None else getattr(orbit, attr)

    row = {
        "time_yr": sim.t,
        "rel_energy_error": abs((sim.energy() - energy_initial) / energy_initial),
        "central_a": value(central_orbit, "a"),
        "central_e": value(central_orbit, "e"),
        "outer_binary_a": value(outer_binary_orbit, "a"),
        "outer_binary_e": value(outer_binary_orbit, "e"),
        "wide_a": value(wide_orbit, "a"),
        "wide_e": value(wide_orbit, "e"),
    }
    row.update(disk_summary)
    return row


def save_snapshot_plot(path, sim, plot_radius, title):
    central_com = sim.com(last=2)
    _, disk_xyz = disk_radii(sim)
    star_xyz = relative_positions(sim, 0, 4, central_com)

    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    if disk_xyz.size:
        ax.scatter(
            disk_xyz[:, 0],
            disk_xyz[:, 1],
            s=0.8,
            c="black",
            alpha=0.35,
            linewidths=0,
        )
    ax.scatter(star_xyz[:2, 0], star_xyz[:2, 1], s=34, c="red", label="central binary")
    ax.scatter(star_xyz[2:, 0], star_xyz[2:, 1], s=34, c="tab:blue", label="outer binary")
    ax.set_xlim(-plot_radius, plot_radius)
    ax.set_ylim(-plot_radius, plot_radius)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x relative to central COM (AU)")
    ax.set_ylabel("y relative to central COM (AU)")
    ax.set_title(title)
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_histogram(path, sim, survivor_radius, title):
    radii, _ = disk_radii(sim)
    radii = radii[np.isfinite(radii)]
    radii = radii[radii <= survivor_radius]

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.hist(radii, bins=80, range=(0.0, survivor_radius), color="black", alpha=0.8)
    ax.axvline(180.0, color="tab:red", linewidth=1.0, label="observed inner edge")
    ax.axvline(260.0, color="tab:blue", linewidth=1.0, label="observed dust ring outer")
    ax.set_xlabel("radius from central COM (AU)")
    ax.set_ylabel("particle count")
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary_csv(path, rows):
    fieldnames = [
        "time_yr",
        "rel_energy_error",
        "central_a",
        "central_e",
        "outer_binary_a",
        "outer_binary_e",
        "wide_a",
        "wide_e",
        "n_within_survivor_radius",
        "inner_edge_p05",
        "outer_edge_p95",
        "median_radius",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def save_time_series(path, rows):
    t = np.array([row["time_yr"] for row in rows], dtype=float)
    energy = np.array([row["rel_energy_error"] for row in rows], dtype=float)
    n_survive = np.array([row["n_within_survivor_radius"] for row in rows], dtype=float)
    inner = np.array([row["inner_edge_p05"] for row in rows], dtype=float)
    outer = np.array([row["outer_edge_p95"] for row in rows], dtype=float)

    fig, axes = plt.subplots(3, 1, figsize=(8.0, 8.0), sharex=True)
    axes[0].semilogy(t, np.where(energy > 0.0, energy, np.nan), "o-", color="black")
    axes[0].set_ylabel("abs(dE/E0)")
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(t, n_survive, "o-", color="black")
    axes[1].set_ylabel("N(r <= limit)")
    axes[1].grid(True, alpha=0.25)

    axes[2].plot(t, inner, "o-", color="tab:red", label="5th percentile")
    axes[2].plot(t, outer, "o-", color="tab:blue", label="95th percentile")
    axes[2].set_ylabel("disk edge proxy (AU)")
    axes[2].set_xlabel("time (yr)")
    axes[2].grid(True, alpha=0.25)
    axes[2].legend(frameon=False)

    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def parse_snapshot_times(args):
    if args.snapshot_times:
        times = [float(t) for t in args.snapshot_times]
    elif args.paper:
        times = [0.0, 44000.0, 1.5e6, 1.5e7]
    else:
        times = [0.0, args.first_periastron_time]
    times = sorted(set(times))
    if times[0] != 0.0:
        times.insert(0, 0.0)
    return times


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run a qualitative Beust section 3.3 GG Tau whfast_hj test."
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--quick", action="store_true", help="Run the early-sculpting case.")
    mode.add_argument("--paper", action="store_true", help="Use Beust Fig. 10 snapshot times.")
    parser.add_argument("--n-disk", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--snapshot-times", nargs="*", type=float)
    parser.add_argument("--plot-radius", type=float, default=1400.0)
    parser.add_argument("--survivor-radius", type=float, default=2000.0)
    parser.add_argument("--dt-divisor", type=float, default=20.0)
    parser.add_argument("--central-a", type=float, default=CENTRAL_A_AU)
    parser.add_argument("--central-e", type=float, default=CENTRAL_E)
    parser.add_argument("--outer-binary-a", type=float, default=OUTER_BINARY_A_AU)
    parser.add_argument("--outer-binary-e", type=float, default=OUTER_BINARY_E)
    parser.add_argument("--wide-periastron", type=float, default=WIDE_PERIASTRON_AU)
    parser.add_argument("--wide-e", type=float, default=WIDE_E)
    parser.add_argument("--wide-i", type=float, default=WIDE_I_DEG)
    parser.add_argument("--outer-binary-rel-i", type=float, default=OUTER_BINARY_REL_I_DEG)
    parser.add_argument("--first-periastron-time", type=float, default=FIRST_WIDE_PERIASTRON_YR)
    parser.add_argument(
        "--allow-large-tree",
        action="store_true",
        help="Allow explicit Python tree generation above the conservative disk limit.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    require_given_tree_api()
    if args.n_disk < 0:
        raise ValueError("--n-disk must be non-negative.")
    if args.n_disk > MAX_EXPLICIT_DISK_PARTICLES and not args.allow_large_tree:
        raise ValueError(
            f"--n-disk={args.n_disk} would build a deeply nested explicit tree. "
            f"Use <= {MAX_EXPLICIT_DISK_PARTICLES} for this prototype, or pass "
            "--allow-large-tree if you intentionally want to test it."
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    snapshot_times = parse_snapshot_times(args)
    sim = setup_ggtau_sim(args)
    tree = make_fixed_tree(sim.N)
    energy_initial = sim.energy()

    print("fixed whfast_hj massive tree: [[1,2],[3,4]]")
    print(f"disk particles: {args.n_disk}")
    print(f"full fixed tree length: {len(tree)} characters")
    print(f"central dt = {sim.dt:.6g} yr")
    print(f"wide a = {wide_semimajor_axis(args.wide_periastron, args.wide_e):.6g} AU")
    print(f"snapshot times: {snapshot_times}")
    print("output dir:", args.output_dir)

    rows = []
    for target_time in snapshot_times:
        if target_time > sim.t:
            sim.integrate(target_time, exact_finish_time=0, given_tree=True, tree=tree)

        row = snapshot_summary(sim, energy_initial, args.survivor_radius)
        rows.append(row)
        tag = f"t{int(round(row['time_yr'])):010d}"
        title = f"GG Tau fixed-tree whfast_hj, t = {row['time_yr']:.0f} yr"
        save_snapshot_plot(
            args.output_dir / f"beust_33_snapshot_{tag}.png",
            sim,
            args.plot_radius,
            title,
        )
        save_histogram(
            args.output_dir / f"beust_33_radial_hist_{tag}.png",
            sim,
            args.survivor_radius,
            title,
        )
        print(
            f"t={row['time_yr']:.3f} yr, "
            f"dE/E={row['rel_energy_error']:.3e}, "
            f"N(r<={args.survivor_radius:g})={row['n_within_survivor_radius']}, "
            f"edge proxies=({row['inner_edge_p05']:.3g}, {row['outer_edge_p95']:.3g}) AU"
        )

    write_summary_csv(args.output_dir / "beust_33_summary.csv", rows)
    save_time_series(args.output_dir / "beust_33_diagnostics.png", rows)
    print("wrote outputs to", args.output_dir)


if __name__ == "__main__":
    main()
