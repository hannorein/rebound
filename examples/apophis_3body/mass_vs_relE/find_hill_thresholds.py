#!/usr/bin/env python3
"""Calibrate mutual-Hill thresholds for the Apophis mass-vs-relE experiment.

This script uses the baseline saved snapshot together with the old
large-mass switch-window sweep to produce candidate constants

    x = start threshold in units of R_H
    a = switch-in threshold in units of R_H
    b = switch-out threshold in units of R_H

The Hill radius is the Earth-Apophis mutual Hill radius around the Sun:

    R_H = ((m_E + m_A) / (3 M_sun))^(1/3) * ((a_E + a_A) / 2)

By default the script picks the best saved `whfast_fake_hj` row
(minimum `max_relE`) from `old_data/large_mass_data/my_sweep.csv`.

The future mass-sweep experiment is intended to use x again as the
outbound stop threshold, so this script does not calibrate a separate y.
"""

from __future__ import annotations

import csv
import os
from ctypes import byref
from pathlib import Path
import sys
import sysconfig
import warnings


SCRIPT_DIR = Path(__file__).resolve().parent
APOPHIS_DIR = SCRIPT_DIR.parent
REPO_ROOT = SCRIPT_DIR.parents[2]

# Make sure this example script imports the local REBOUND package.
sys.path.insert(0, str(REPO_ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(APOPHIS_DIR / ".mplconfig"))

# Match REBOUND's expected shared-library location when running from source.
suffix = sysconfig.get_config_var("EXT_SUFFIX") or ".so"
expected_lib = REPO_ROOT / f"librebound{suffix}"
source_lib = REPO_ROOT / "src" / "librebound.so"
if not expected_lib.exists() and source_lib.exists():
    expected_lib.symlink_to(source_lib.relative_to(REPO_ROOT))

import rebound


GAUSSIAN_K = 0.01720209895
STANDARD_ORDER = ("sun", "earth", "asteroid")
FAKE_HJ_ORDER = ("earth", "asteroid", "sun")

# One-off calibration inputs for the baseline Apophis setup.
# Change these here if you ever want to calibrate from a different baseline.
SNAPSHOT = APOPHIS_DIR / "find_close_encounter" / "apophis_1yr_before.bin"
SUMMARY_CSV = APOPHIS_DIR / "old_data" / "large_mass_data" / "my_sweep.csv"
MODE = "whfast_fake_hj"
CASE = None  # Set to an integer like 3 to force a specific saved sweep case.
DURATION_DAYS = 365.25 * 2.0
WHFAST_DT = 10.0


def load_summary_row(summary_csv: Path, mode: str, case: int | None) -> dict[str, str]:
    # Read the saved sweep table and keep only rows for the requested mode.
    rows = list(csv.DictReader(summary_csv.open()))
    mode_rows = [row for row in rows if row["mode"] == mode]
    if not mode_rows:
        raise ValueError(f"No rows found for mode {mode!r} in {summary_csv}.")

    if case is not None:
        # If a specific case is requested, use exactly that saved row.
        matches = [row for row in mode_rows if int(row["case"]) == case]
        if not matches:
            raise ValueError(f"No row found for mode {mode!r} and case={case} in {summary_csv}.")
        return matches[0]

    # Otherwise use the saved row with the smallest maximum relative energy error.
    return min(mode_rows, key=lambda row: float(row["max_relE"]))


def load_sim(snapshot: Path) -> rebound.Simulation:
    # Load the binary snapshot, then restore REBOUND's function pointers.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="You have to reset function pointers after creating a reb_simulation struct with a binary file.",
            category=RuntimeWarning,
        )
        sim = rebound.Simulation(str(snapshot))
    rebound.clibrebound.reb_simulation_reset_function_pointers(byref(sim))

    # Restore the expected hashes so we can refer to particles by name.
    for particle, key in zip(sim.particles[:3], STANDARD_ORDER):
        particle.hash = key
    return sim


def configure_whfast(sim: rebound.Simulation, whfast_dt: float) -> None:
    # Rebuild the same WHFast setup used in the old baseline comparison.
    sim.N_active = 3
    sim.force_is_velocity_dependent = 0
    sim.G = GAUSSIAN_K * GAUSSIAN_K
    sim.reset_integrator()
    sim.integrator = "whfast"
    sim.dt = whfast_dt
    sim.ri_whfast.safe_mode = 1
    sim.ri_whfast.corrector = 11
    sim.ri_whfast.coordinates = "jacobi"
    sim.ri_whfast.jacobi_ordered_warning = -1


def clone_with_order(source_sim: rebound.Simulation, order: tuple[str, ...], whfast_dt: float) -> rebound.Simulation:
    # Rebuild a simulation from the current inertial state but with a different
    # particle ordering. This is how the helper mimics fake HJ.
    new_sim = rebound.Simulation()
    new_sim.t = source_sim.t
    new_sim.dt = source_sim.dt
    new_sim.softening = source_sim.softening
    new_sim.testparticle_type = source_sim.testparticle_type
    for key in order:
        new_sim.add(source_sim.particles[key].copy())
        new_sim.particles[-1].hash = key
    configure_whfast(new_sim, whfast_dt)
    return new_sim


def advance_whfast_by_delta(sim: rebound.Simulation, delta_days: float, whfast_dt: float) -> None:
    # Advance with repeated WHFast steps so the replay follows the same
    # fixed-step structure as the old baseline branch.
    if delta_days < -1e-14:
        raise ValueError(f"Cannot step backward by {delta_days} days.")
    if delta_days <= 1e-14:
        return

    full_steps = int(delta_days // whfast_dt)
    remainder = delta_days - whfast_dt * full_steps

    for _ in range(full_steps):
        sim.dt = whfast_dt
        sim.step()

    if remainder > 1e-14:
        sim.dt = remainder
        sim.step()


def earth_asteroid_distance(sim: rebound.Simulation) -> float:
    # Physical encounter distance used for the future Hill-based trigger.
    particles = sim.particles
    dx = particles["asteroid"].x - particles["earth"].x
    dy = particles["asteroid"].y - particles["earth"].y
    dz = particles["asteroid"].z - particles["earth"].z
    return (dx * dx + dy * dy + dz * dz) ** 0.5


def mutual_hill_radius(sim: rebound.Simulation) -> tuple[float, float, float]:
    # Compute the chosen mutual Hill radius for the Earth-Apophis pair.
    # Returning a_E and a_A as well makes the printed diagnostics easier to read.
    particles = sim.particles
    sun = particles["sun"]
    earth = particles["earth"]
    asteroid = particles["asteroid"]
    earth_orbit = earth.orbit(primary=sun)
    asteroid_orbit = asteroid.orbit(primary=sun)
    rh = ((earth.m + asteroid.m) / (3.0 * sun.m)) ** (1.0 / 3.0) * ((earth_orbit.a + asteroid_orbit.a) / 2.0)
    return rh, earth_orbit.a, asteroid_orbit.a


def print_event(name: str, sim: rebound.Simulation) -> None:
    # Convenience printer for one state along the reconstructed baseline branch.
    distance = earth_asteroid_distance(sim)
    rh, a_earth, a_asteroid = mutual_hill_radius(sim)
    print(f"{name}:")
    print(f"  t [day]      = {sim.t:.16f}")
    print(f"  d_EA [AU]    = {distance:.16e}")
    print(f"  R_H [AU]     = {rh:.16e}")
    print(f"  d_EA / R_H   = {distance / rh:.16e}")
    print(f"  a_E [AU]     = {a_earth:.16e}")
    print(f"  a_A [AU]     = {a_asteroid:.16e}")


def main() -> None:
    # Step 1: locate the fixed baseline inputs.
    snapshot = SNAPSHOT.resolve()
    summary_csv = SUMMARY_CSV.resolve()
    if not snapshot.exists():
        raise FileNotFoundError(f"Snapshot file not found: {snapshot}")
    if not summary_csv.exists():
        raise FileNotFoundError(f"Sweep summary CSV not found: {summary_csv}")

    # Step 2: choose the saved sweep row that defines the baseline switch window.
    row = load_summary_row(summary_csv, MODE, CASE)
    switch_in_time = float(row["actual_switch_start_days"])
    switch_out_time = float(row["actual_switch_end_days"])

    # Step 3: load the snapshot and configure it like the old WHFast run.
    sim = load_sim(snapshot)
    configure_whfast(sim, WHFAST_DT)
    start_time = sim.t

    # Step 4: x is the old simulation start distance in Hill units.
    # In the future mass sweep, this same x will be reused as the outbound stop
    # threshold as well.
    # x: old restart state before any switching.
    x_distance = earth_asteroid_distance(sim)
    x_rh, _, _ = mutual_hill_radius(sim)
    x_value = x_distance / x_rh

    # Step 5: advance to the saved switch-in time and measure a on the standard
    # WHFast branch, just before switching into fake HJ.
    # a: state immediately before switching into fake-HJ ordering.
    advance_whfast_by_delta(sim, switch_in_time - start_time, WHFAST_DT)
    a_distance = earth_asteroid_distance(sim)
    a_rh, _, _ = mutual_hill_radius(sim)
    a_value = a_distance / a_rh

    # Step 6: reorder into the fake-HJ layout, continue to the saved switch-out
    # time, and measure b on that switched branch.
    # b: follow the actual fake-HJ branch until the saved switch-out time.
    sim = clone_with_order(sim, FAKE_HJ_ORDER, WHFAST_DT)
    advance_whfast_by_delta(sim, switch_out_time - switch_in_time, WHFAST_DT)
    b_distance = earth_asteroid_distance(sim)
    b_rh, _, _ = mutual_hill_radius(sim)
    b_value = b_distance / b_rh

    # Step 7: emit only the constants to paste into the main mass-sweep script.
    print(f"X_HILL = {x_value:.16f}")
    print(f"A_HILL = {a_value:.16f}")
    print(f"B_HILL = {b_value:.16f}")


if __name__ == "__main__":
    main()
