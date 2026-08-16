#!/usr/bin/env python3
"""Compare a retained WHFast/HJ case to an adaptive IAS15 trajectory."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from search_fixed_hierarchy import DEFAULT_OUTPUT, build_simulation

import matplotlib.pyplot as plt


def state(sim):
    return np.asarray(
        [[p.x, p.y, p.z, p.vx, p.vy, p.vz] for p in sim.particles], dtype=float
    )


def hierarchy_vectors(states, masses):
    left_mass = masses[0] + masses[1]
    right_mass = masses[2] + masses[3]
    left_com = (masses[0] * states[:, 0] + masses[1] * states[:, 1]) / left_mass
    right_com = (masses[2] * states[:, 2] + masses[3] * states[:, 3]) / right_mass
    return np.stack(
        (
            states[:, 1] - states[:, 0],
            states[:, 3] - states[:, 2],
            right_com - left_com,
        ),
        axis=1,
    )


def normalized_rms(candidate, reference, columns):
    difference = candidate[..., columns] - reference[..., columns]
    scale = np.sqrt(np.mean(reference[..., columns] ** 2))
    return float(np.sqrt(np.mean(difference**2)) / scale)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_OUTPUT / "rank_01")
    parser.add_argument("--samples", type=int, default=1001)
    parser.add_argument("--ias15-epsilon", type=float, default=1.0e-12)
    parser.add_argument("--duration-multiplier", type=float, default=1.0)
    return parser.parse_args()


def main():
    args = parse_args()
    with (args.input / "parameters.json").open() as stream:
        config = json.load(stream)
    if args.duration_multiplier <= 0.0:
        raise ValueError("duration multiplier must be positive")
    config["duration"] *= args.duration_multiplier
    config["n_steps"] = int(math.ceil(config["duration"] / config["dt"]))

    sims = {
        "whfast": build_simulation(config, "whfast"),
        "whfast_hj": build_simulation(config, "whfast_hj"),
        "ias15": build_simulation(config, "ias15"),
    }
    sims["ias15"].integrator.epsilon = args.ias15_epsilon
    initial_energy = {mode: sim.energy() for mode, sim in sims.items()}
    sample_steps = np.unique(
        np.linspace(0, config["n_steps"], args.samples, dtype=np.int64)
    )
    histories = {mode: [] for mode in sims}
    energy_errors = {mode: [] for mode in sims}
    times = []
    previous = 0
    for step in sample_steps:
        delta = int(step - previous)
        if delta:
            sims["whfast"].steps(delta)
            sims["whfast_hj"].steps(delta)
            target = step * config["dt"]
            sims["ias15"].integrate(target, exact_finish_time=1)
        times.append(step * config["dt"])
        for mode, sim in sims.items():
            histories[mode].append(state(sim))
            energy_errors[mode].append(
                max(
                    abs((sim.energy() - initial_energy[mode]) / initial_energy[mode]),
                    np.finfo(float).eps,
                )
            )
        previous = step

    masses = np.asarray(config["masses"])
    vectors = {
        mode: hierarchy_vectors(np.asarray(history), masses)
        for mode, history in histories.items()
    }
    reference = vectors["ias15"]
    report = {
        "id": config["id"],
        "ias15_epsilon": args.ias15_epsilon,
        "samples": len(sample_steps),
        "duration_multiplier": args.duration_multiplier,
        "ias15_max_energy_error": float(np.max(energy_errors["ias15"])),
    }
    for mode in ("whfast", "whfast_hj"):
        report[mode] = {
            "hierarchy_position_rms_error": normalized_rms(vectors[mode], reference, slice(0, 3)),
            "hierarchy_velocity_rms_error": normalized_rms(vectors[mode], reference, slice(3, 6)),
            "sampled_max_energy_error": float(np.max(energy_errors[mode])),
        }
    with (args.input / "ias15_validation.json").open("w") as stream:
        json.dump(report, stream, indent=2)

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    colors = {"whfast": "#d95f02", "whfast_hj": "#1b9e77", "ias15": "#3b5bdb"}
    labels = {"whfast": "WHFast", "whfast_hj": "fixed HJ", "ias15": "IAS15"}
    normalized_time = np.asarray(times) / config["reference_period"]
    for mode in ("whfast", "whfast_hj", "ias15"):
        ax.semilogy(normalized_time, energy_errors[mode], color=colors[mode], lw=1, label=labels[mode])
    ax.set_xlabel(r"time / $P_{\rm ref}$")
    ax.set_ylabel(r"$|E(t)-E(0)|/|E(0)|$")
    ax.set_title(f"IAS15 audit: {config['id']}")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.input / "ias15_validation.png", dpi=180)
    plt.close(fig)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
