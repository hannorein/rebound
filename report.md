# REBOUND Project Deep-Dive Report

## 1. What this project is
REBOUND is a mature N-body framework (C core + Python bindings) for astrophysical dynamics. It provides multiple integrators (symplectic, hybrid, adaptive), collision handling, simulation archives, and a large test/example suite.

Your research goal:
- Use **hierarchical Jacobi coordinates** in a **hybrid symplectic integrator** for planetary systems with close encounters.
- Part 1: implement algorithm.
- Part 2: validate accuracy in astrophysical scenarios.
- Part 3: optimize and benchmark vs Mercury/Samba.

The repository already contains almost all infrastructure needed for this program, but not your specific coordinate/splitting upgrade in hybrid mode.

## 2. High-level architecture

### 2.1 Layers
- **Numerical core (C):** `src/*.c`, public API in `src/rebound.h`.
- **Python wrapper (ctypes):** `rebound/*.py`, integrator-specific wrappers in `rebound/integrators/*.py`.
- **Validation and examples:** `rebound/tests/*.py`, `examples/*`, `python_examples/*`, `ipython_examples/*`.
- **Docs:** `docs/*.md`.

### 2.2 Main runtime flow
Single timestep entry point is `reb_simulation_step()` in `src/rebound.c`:
1. Optional pre-step modification callbacks.
2. `reb_integrator_part1(r)` dispatch in `src/integrator.c`.
3. Tree/boundary updates if needed.
4. Gravity + additional forces (`reb_calculate_acceleration`, callbacks).
5. `reb_integrator_part2(r)` dispatch.
6. Optional post-step modification callbacks.
7. Collision search (`reb_collision_search`).
8. Wallclock bookkeeping and `steps_done` increment.

This split (`part1`/`part2`) is the core integration contract for all built-in integrators.

## 3. Integrator system and where close-encounter logic lives

### 3.1 Dispatch and reset
- Dispatch: `src/integrator.c` (`reb_integrator_part1/part2`, synchronize, reset).
- Integrator selection enum in `src/rebound.h` (`REB_INTEGRATOR_*`).

### 3.2 WHFast (symplectic baseline)
- Implementation: `src/integrator_whfast.c`.
- Key features already present:
  - coordinate choices: Jacobi, democratic heliocentric, WHDS, barycentric.
  - Kepler/interactions/jump/com operators.
  - kernels + symplectic correctors.
  - Jacobi hierarchy validation hooks via orbit hierarchy.
- Important checks in `reb_integrator_whfast_init()`:
  - non-default kernels require Jacobi.
  - variational eqs only with Jacobi + default kernel.
  - warnings if particles are not Jacobi-ordered or not central-body hierarchy.

### 3.3 MERCURIUS (hybrid symplectic)
- Implementation: `src/integrator_mercurius.c`.
- Current strategy:
  - Normal evolution: WHFast-like split in **democratic heliocentric internals**.
  - Encounter handling: predict encounters (`reb_mercurius_encounter_predict`), then integrate encounter subset with IAS15 (`reb_mercurius_encounter_step`).
  - Uses switching function `L(d,dcrit)` and `dcrit` criteria including Hill/velocity/radius terms.
- Internal coordinate transforms:
  - `reb_integrator_mercurius_inertial_to_dh`
  - `reb_integrator_mercurius_dh_to_inertial`
- Gravity mode set to `REB_GRAVITY_MERCURIUS`.

### 3.4 TRACE (hybrid reversible)
- Implementation: `src/integrator_trace.c`.
- Current strategy:
  - Uses WHFast operators for long-term regions.
  - Detects close encounters/pericenter conditions (`S`, `S_peri`).
  - Time-reversible step rejection/redo if new encounters appear post-step.
  - Encounter integration via BS or IAS15 depending on `peri_mode`.
- Internal coordinate transforms:
  - `reb_integrator_trace_inertial_to_dh`
  - `reb_integrator_trace_dh_to_inertial`
- Gravity mode set to `REB_GRAVITY_TRACE`.

## 4. Coordinate infrastructure status vs your research goal

### 4.1 Existing transformation library
`src/transformations.c` already implements conversions for:
- inertial <-> Jacobi
- inertial <-> democratic heliocentric
- inertial <-> WHDS
- inertial <-> barycentric

No hierarchical-Jacobi hybrid transform path exists yet for MERCURIUS/TRACE.

### 4.2 Orbit hierarchy support already exists
- Core tree builder and checks: `src/orbit_hierarchy.c`.
- Public API in `src/rebound.h`:
  - `reb_orbit_hierarchy_create_from_simulation`
  - `reb_orbit_hierarchy_is_jacobi`
  - `reb_orbit_hierarchy_is_jacobi_ordered`
- Python wrapper: `rebound/orbit_hierarchy.py`.
- Tests: `rebound/tests/test_orbit_hierarchy.py`.

This is a strong foundation for hierarchical coordinate work.

### 4.3 Current hybrid limitation (critical finding)
MERCURIUS and TRACE currently operate with democratic-heliocentric style internals and encounter switching, not hierarchical Jacobi coordinate evolution. Your target method is therefore a genuine extension, not a small parameter change.

## 5. Data model essentials
- Main struct: `struct reb_simulation` in `src/rebound.h`.
- Integrator state structs embedded in simulation:
  - `ri_whfast`, `ri_mercurius`, `ri_trace`, etc.
- Function pointer extension points:
  - additional forces, heartbeat, collision resolve, pre/post timestep hooks.

Any new hybrid integrator mode should either:
1. extend `ri_mercurius`/`ri_trace` with coordinate/splitting mode fields, or
2. introduce a new integrator enum + struct if behavior diverges significantly.

## 6. Testing and CI maturity
- Python CI: `.github/workflows/python.yml` (unittest over `rebound/tests`, example scripts).
- C CI: `.github/workflows/c.yml` (build shared library + compile all C examples).
- Existing tests cover WHFast, MERCURIUS, TRACE (including simulation archive/restart, encounter conditions, energy behavior, ordering invariance).

This suite is sufficient to support your Part 2 accuracy program with additional scenario-specific tests.

## 7. Concrete implementation map for your Part 1

### 7.1 Most likely code entry points
- `src/integrator_mercurius.c`
- `src/integrator_trace.c`
- `src/transformations.c`
- `src/gravity.c` (if new splitting requires modified force partition)
- `src/rebound.h` (new mode flags/state)
- `rebound/integrators/mercurius.py` and/or `trace.py` (expose new knobs in Python)
- `docs/integrators.md` (document behavior)

### 7.2 Recommended strategy
- Add a coordinate/splitting mode flag to one hybrid integrator first (MERCURIUS is usually simpler than TRACE’s reversibility constraints).
- Implement inertial <-> hierarchical-Jacobi transform path (can leverage orbit hierarchy tree for ordering).
- Implement Hamiltonian split terms for new mode while preserving old defaults.
- Keep backward-compatible defaults to avoid breaking existing test corpus.

### 7.3 Invariants to preserve
- Symplectic segment should remain time-centered correctly around encounter switching.
- Synchronization semantics (`is_synchronized`, recalc flags) must stay consistent.
- Collision map and encounter map indices must stay valid across coordinate conversions.
- Simulation archive compatibility: if you add persistent struct fields, ensure binary field descriptor support (`src/output.c`/serialization paths).

## 8. Part 2 (accuracy testing) blueprint
Use existing test and example patterns to add:
- Long-term secular drift tests (energy/angular momentum) for non-encounter systems.
- Controlled close-encounter tests with known IAS15 reference runs.
- Multi-planet packed systems (already similar fixtures in `test_mercurius.py`, `test_trace*.py`).
- Hierarchical systems with moons/circumbinary-style architecture to validate hierarchy handling.

Metrics:
- Relative energy drift over fixed physical time.
- Encounter trajectory mismatch vs high-accuracy reference.
- Reversibility checks (especially if extending TRACE-like approach).

## 9. Part 3 (performance/benchmarking) blueprint
- Baseline timing harnesses can be cloned from existing timing-oriented tests in trace/mercurius suites.
- Compare:
  - existing MERCURIUS/TRACE mode
  - new hierarchical-Jacobi hybrid mode
  - IAS15 reference (accuracy baseline)
- For external Mercury/Samba comparison, add reproducible initial-condition export/import scripts and standardized runtime/accuracy metrics.

## 10. Architecture confidence statement
I am confident in the architecture assessment:
- The simulation stepping contract and dispatch structure are clear and stable.
- Hybrid close-encounter behavior is localized and explicit in `integrator_mercurius.c` and `integrator_trace.c`.
- Coordinate transform capabilities are comprehensive but currently do not include a hierarchical-Jacobi hybrid path.
- Orbit hierarchy utilities are present and reusable for your target design.

This codebase is well-positioned for your three-phase research plan with incremental, test-driven extension.
