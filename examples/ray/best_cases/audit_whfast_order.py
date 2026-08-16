#!/usr/bin/env python3
"""Audit all 24 ordinary-WHFast particle orders for one retained case."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from pathlib import Path

import numpy as np

from search_fixed_hierarchy import DEFAULT_OUTPUT, build_simulation, evaluate

import rebound


NAMES = ("A1", "A2", "B1", "B2")


def reordered_whfast(config: dict, order: tuple[str, ...]):
    original = build_simulation(config, "whfast")
    states = {}
    for name, particle in zip(NAMES, original.particles):
        states[name] = (
            particle.m,
            particle.x,
            particle.y,
            particle.z,
            particle.vx,
            particle.vy,
            particle.vz,
        )
    sim = rebound.Simulation()
    sim.G = 1.0
    for name in order:
        m, x, y, z, vx, vy, vz = states[name]
        sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, name=name)
    sim.integrator = "whfast"
    sim.dt = config["dt"]
    return sim


def maximum_energy_error(sim, n_steps: int):
    initial = sim.energy()
    maximum = np.finfo(float).eps
    for _ in range(n_steps):
        sim.steps(1)
        error = abs((sim.energy() - initial) / initial)
        if not np.isfinite(error):
            return np.nan
        maximum = max(maximum, float(error))
    return maximum


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_OUTPUT / "rank_04")
    return parser.parse_args()


def main():
    args = parse_args()
    with (args.input / "parameters.json").open() as stream:
        config = json.load(stream)

    hj_error = evaluate(config, samples=2)["whfast_hj_max_error"]
    rows = []
    for order in itertools.permutations(NAMES):
        error = maximum_energy_error(reordered_whfast(config, order), config["n_steps"])
        rows.append(
            {
                "particle_order": "-".join(order),
                "whfast_max_error": error,
                "fixed_hj_max_error": hj_error,
                "improvement": error / hj_error,
            }
        )
    rows.sort(key=lambda row: row["whfast_max_error"])
    with (args.input / "whfast_order_audit.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print("Best ordinary-WHFast order:", rows[0])
    print("Worst ordinary-WHFast order:", rows[-1])
    print("Median improvement:", np.median([row["improvement"] for row in rows]))


if __name__ == "__main__":
    main()
