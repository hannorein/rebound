# Hill-radius definition:
#     R_H = ((m_E + m_A) / (3 M_sun))^(1/3) * ((a_E + a_A) / 2)
#     m_E, m_A = Earth and Apophis masses
#     M_sun    = Sun mass
#     a_E, a_A = Earth and Apophis semimajor axes
#
# Energy-error definition:
#   relE(t) = (E(t) - E0) / E0
# where
#   E(t) = total simulation energy at time t
#   E0   = total simulation energy at the common restart time t_start
# Quantity to plot:
#   max_t |relE(t)|
# measured over the integration interval [t_start, t_stop].

from __future__ import annotations

import sys
import sysconfig
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
APOPHIS_DIR = SCRIPT_DIR.parent
REPO_ROOT = SCRIPT_DIR.parents[2]

sys.path.insert(0, str(REPO_ROOT))

suffix = sysconfig.get_config_var("EXT_SUFFIX") or ".so"
expected_lib = REPO_ROOT / f"librebound{suffix}"
source_lib = REPO_ROOT / "src" / "librebound.so"
if not expected_lib.exists() and source_lib.exists():
    expected_lib.symlink_to(source_lib.relative_to(REPO_ROOT))

import rebound

OUTPUT_DIR = SCRIPT_DIR / "outputs"

GAUSSIAN_K = 0.01720209895
IAS15_INITIAL_DT_DAYS = 1.0

SUN_MASS = 1.0
EARTH_MASS_0 = 3.0034896149156e-6
APOPHIS_MASS_0 = 3.0034896149156e-6

BASE_POSITIONS = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.9224, 0.0, 0.0),
)

BASE_VELOCITIES = (
    (0.0, 0.0, 0.0),
    (0.0, GAUSSIAN_K, 0.0),
    (0.0, 0.0185, 0.0),
)

X_HILL = 18.0
A_HILL = 3.6
B_HILL = 4.38
LAMBDA_VALUES = (1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0)
STANDARD_ORDER = ("sun", "earth", "asteroid")
FAKE_HJ_ORDER = ("earth", "asteroid", "sun")

def make_base_simulation() -> rebound.Simulation:
    """Build the base IAS15 simulation from the original initial conditions."""
    sim = rebound.Simulation()
    sim.G = GAUSSIAN_K * GAUSSIAN_K
    sim.dt = IAS15_INITIAL_DT_DAYS
    sim.force_is_velocity_dependent = 0
    sim.integrator = "ias15"

    masses = (SUN_MASS, EARTH_MASS_0, APOPHIS_MASS_0)
    for key, pos, vel, mass in zip(STANDARD_ORDER, BASE_POSITIONS, BASE_VELOCITIES, masses):
        sim.add(
            m=mass,
            x=pos[0],
            y=pos[1],
            z=pos[2],
            vx=vel[0],
            vy=vel[1],
            vz=vel[2],
        )

    sim.move_to_com()
    sim.N_active = 3
    for particle, key in zip(sim.particles[:3], STANDARD_ORDER):
        particle.hash = key
    return sim


def apply_mass_scale(
    sim: rebound.Simulation,
    lambda_value: float,
) -> None:
    """Apply step 1 to a simulation branch using the baseline masses."""
    sim.particles["earth"].m = EARTH_MASS_0 * lambda_value
    sim.particles["asteroid"].m = APOPHIS_MASS_0 * lambda_value


def earth_asteroid_distance(sim: rebound.Simulation) -> float:
    """Return the Earth-Apophis separation in AU."""
    particles = sim.particles
    dx = particles["asteroid"].x - particles["earth"].x
    dy = particles["asteroid"].y - particles["earth"].y
    dz = particles["asteroid"].z - particles["earth"].z
    return (dx * dx + dy * dy + dz * dz) ** 0.5


def mutual_hill_radius(sim: rebound.Simulation) -> float:
    """Return the Earth-Apophis mutual Hill radius around the Sun."""
    particles = sim.particles
    sun = particles["sun"]
    earth = particles["earth"]
    asteroid = particles["asteroid"]
    earth_orbit = earth.orbit(primary=sun)
    asteroid_orbit = asteroid.orbit(primary=sun)
    return ((earth.m + asteroid.m) / (3.0 * sun.m)) ** (1.0 / 3.0) * ((earth_orbit.a + asteroid_orbit.a) / 2.0)


def configure_whfast(sim: rebound.Simulation) -> None:
    """Configure a simulation branch for the WHFast fake-HJ experiment."""
    sim.N_active = 3
    sim.force_is_velocity_dependent = 0
    sim.G = GAUSSIAN_K * GAUSSIAN_K
    sim.reset_integrator()
    sim.integrator = "whfast"
    sim.ri_whfast.safe_mode = 1
    sim.ri_whfast.corrector = 11
    sim.ri_whfast.coordinates = "jacobi"
    sim.ri_whfast.jacobi_ordered_warning = -1


def clone_with_order(
    source_sim: rebound.Simulation,
    order: tuple[str, ...],
) -> rebound.Simulation:
    """Rebuild the current inertial state with a different particle order."""
    new_sim = rebound.Simulation()
    new_sim.t = source_sim.t
    new_sim.dt = source_sim.dt
    new_sim.softening = source_sim.softening
    new_sim.testparticle_type = source_sim.testparticle_type
    for key in order:
        new_sim.add(source_sim.particles[key].copy())
        new_sim.particles[-1].hash = key
    configure_whfast(new_sim)
    return new_sim


def run_ias15_until_x_hill_exit(
    sim: rebound.Simulation,
) -> tuple[float, rebound.Simulation, float | None, float | None, float, float]:
    """Run IAS15 through one close encounter and return a step-end restart snapshot plus absolute times."""
    previous_time = sim.t
    previous_distance = earth_asteroid_distance(sim)
    previous_rh = mutual_hill_radius(sim)
    previous_x = previous_distance - X_HILL * previous_rh
    previous_a = previous_distance - A_HILL * previous_rh
    previous_b = previous_distance - B_HILL * previous_rh
    sim_start = None
    t_switch_in = None
    closest_time = None
    t_switch_out = None
    t_stop = None
    max_dt = 0.0
    relE0 = None
    max_relE = 0.0
    approaching_after_start = False

    def reset_relE_tracker() -> None:
        nonlocal relE0, max_relE
        relE0 = sim.energy()
        max_relE = 0.0

    while True:
        sim.step()
        current_distance = earth_asteroid_distance(sim)
        current_rh = mutual_hill_radius(sim)
        current_x = current_distance - X_HILL * current_rh
        current_a = current_distance - A_HILL * current_rh
        current_b = current_distance - B_HILL * current_rh
        max_dt = max(max_dt, sim.dt_last_done)

        if closest_time is None:
            if previous_x > 0.0 and current_x <= 0.0:
                sim_start = sim.copy()
                t_switch_in = sim.t if current_a <= 0.0 else None
                approaching_after_start = current_distance < previous_distance
                reset_relE_tracker()

            if sim_start is not None and t_switch_in is None and previous_a > 0.0 and current_a <= 0.0:
                t_switch_in = sim.t

            if sim_start is not None:
                if current_distance < previous_distance:
                    approaching_after_start = True
                elif approaching_after_start and current_distance > previous_distance:
                    closest_time = previous_time

        if closest_time is not None:
            if sim_start is not None and t_switch_out is None and previous_b < 0.0 and current_b >= 0.0:
                t_switch_out = sim.t

            if previous_x < 0.0 and current_x >= 0.0:
                if relE0 is None or sim_start is None:
                    raise RuntimeError("Found outbound X_HILL crossing before establishing sim_start.")
                t_stop = sim.t
                relE = (sim.energy() - relE0) / relE0
                max_relE = max(max_relE, abs(relE))
                return max_dt, sim_start, t_switch_in, t_switch_out, t_stop, max_relE

        if relE0 is not None and sim_start is not None:
            relE = (sim.energy() - relE0) / relE0
            max_relE = max(max_relE, abs(relE))

        previous_time = sim.t
        previous_distance = current_distance
        previous_x = current_x
        previous_a = current_a
        previous_b = current_b


def run_WHfast_fake_HJ(
    sim_start: rebound.Simulation,
    t_switch_in: float | None,
    t_switch_out: float | None,
    t_stop: float,
) -> float:
    """Return max |relE| for the WHFast + fake-HJ branch started from sim_start."""
    sim = sim_start
    if sim.dt <= 0.0:
        raise ValueError("sim_start.dt must be set to the WHFast timestep before calling run_WHfast_fake_HJ().")

    configure_whfast(sim)
    e0 = sim.energy()
    max_relE = 0.0
    current_order = STANDARD_ORDER

    while True:
        if current_order == STANDARD_ORDER and t_switch_in is not None and sim.t >= t_switch_in:
            sim = clone_with_order(sim, FAKE_HJ_ORDER)
            current_order = FAKE_HJ_ORDER

        if current_order == FAKE_HJ_ORDER and t_switch_out is not None and sim.t >= t_switch_out:
            sim = clone_with_order(sim, STANDARD_ORDER)
            current_order = STANDARD_ORDER

        if sim.t >= t_stop:
            break

        sim.step()
        relE = (sim.energy() - e0) / e0
        max_relE = max(max_relE, abs(relE))

    return max_relE


def make_plot(results: list[tuple[float, float, float]]) -> Path:
    """Plot max |relE| versus lambda."""
    import matplotlib.pyplot as plt

    output_path = OUTPUT_DIR / "mass_vs_relE.png"
    lambdas = [row[0] for row in results]
    relE_ias15 = [row[1] for row in results]
    relE_wh_hj = [row[2] for row in results]

    fig, ax = plt.subplots()
    ax.loglog(lambdas, relE_ias15, "o-", label="IAS15")
    ax.loglog(lambdas, relE_wh_hj, "s-", label="WHFast + fake HJ")
    ax.set_xlabel("lambda")
    ax.set_ylabel("max |relE|")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return output_path


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results = []

    for lambda_value in LAMBDA_VALUES:
        sim = make_base_simulation()
        apply_mass_scale(sim, lambda_value)
        max_dt, sim_start, t_switch_in, t_switch_out, t_stop, max_relE_IAS15 = run_ias15_until_x_hill_exit(sim)
        sim_start.dt = max_dt
        max_relE_WH_HJ = run_WHfast_fake_HJ(sim_start, t_switch_in, t_switch_out, t_stop)
        results.append((lambda_value, max_relE_IAS15, max_relE_WH_HJ))

    plot_path = make_plot(results)
    print(f"Wrote plot to {plot_path}")


if __name__ == "__main__":
    main()
