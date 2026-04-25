# AGENTS.md

## Purpose
This repository is REBOUND (C core + Python bindings) for astrophysical N-body simulations. Current local research focus is:
- hierarchical Jacobi coordinates in hybrid symplectic close-encounter integrators.

## Core Architecture (must know)
- Main C API and structs: `src/rebound.h`
- Main simulation loop: `src/rebound.c` (`reb_simulation_step`)
- Integrator dispatch/reset/sync: `src/integrator.c`
- Coordinate transforms: `src/transformations.c`
- Gravity modes: `src/gravity.c`
- Collision search: `src/collision.c`
- Orbit hierarchy utilities: `src/orbit_hierarchy.c`, `rebound/orbit_hierarchy.py`

## Integrator-specific files
- WHFast: `src/integrator_whfast.c`, wrapper `rebound/integrators/whfast.py`
- MERCURIUS: `src/integrator_mercurius.c`, wrapper `rebound/integrators/mercurius.py`
- TRACE: `src/integrator_trace.c`, wrapper `rebound/integrators/trace.py`

## Current coordinate reality
- WHFast/SABA support Jacobi and related coordinates.
- MERCURIUS/TRACE currently use democratic-heliocentric style internals + encounter switching.
- No hierarchical-Jacobi hybrid mode exists yet.

## If implementing hierarchical-Jacobi hybrid mode
1. Start in `src/integrator_mercurius.c` (simpler than TRACE reversibility path).
2. Add mode/state fields in `src/rebound.h` integrator struct.
3. Add transform support in `src/transformations.c` (or dedicated helper file if cleaner).
4. Reuse orbit hierarchy tree (`reb_orbit_hierarchy_*`) to enforce/validate hierarchy ordering.
5. Keep default behavior backward-compatible.
6. Update Python wrapper to expose new settings.
7. Update docs in `docs/integrators.md`.
8. Add tests in `rebound/tests/` mirroring style of `test_mercurius.py`, `test_trace*.py`, `test_orbit_hierarchy.py`.

## Testing and CI
- Python CI runs `python -m unittest discover -s rebound/tests/ -v`.
- C CI compiles all C examples under `examples/`.
- Before major changes, at least run focused tests:
  - `rebound/tests/test_orbit_hierarchy.py`
  - `rebound/tests/test_whfast.py`
  - `rebound/tests/test_mercurius.py`
  - `rebound/tests/test_trace*.py`

## Engineering constraints
- Preserve integrator synchronization invariants (`is_synchronized`, recalc flags).
- Preserve encounter map/index correctness under coordinate changes.
- Preserve simulation archive compatibility when adding persistent fields.
- Do not change default numerical behavior unless explicitly intended and tested.

## Useful docs/examples
- Integrator docs: `docs/integrators.md`
- Reference frame docs: `docs/simulationreferenceframes.md`
- Hybrid C example: `examples/closeencounter_hybrid/problem.c`
- TRACE notebook: `ipython_examples/HybridIntegrationsWithTRACE.ipynb`
