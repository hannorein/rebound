# Fixed-hierarchy HJ energy-error search

This experiment searches stable 2+2 quadruples for cases where the maximum
relative energy error of fixed-tree WHFast-HJ is at least 10,000 times smaller
than ordinary WHFast's error. Both integrations use identical Cartesian initial
conditions, timestep, duration, and output cadence.

The fixed hierarchy is `[[A1,A2],[B1,B2]]`, supplied to the integrator as
`[[1,2],[3,4]]`. This topology is intentionally important: ordinary WHFast's
linear Jacobi chain cannot simultaneously represent both tight binaries as
independent Kepler problems.

Run a short search:

```bash
python examples/ray/hj_energy_search/search_fixed_hierarchy.py \
  --candidates 40 --stop-after-target
```

Run a broader search:

```bash
python examples/ray/hj_energy_search/search_fixed_hierarchy.py \
  --candidates 160
```

Validate the retained five over a timestep sequence:

```bash
python examples/ray/hj_energy_search/validate_best.py
```

Audit all 24 ordinary-WHFast particle orderings for the symmetric retained
quadruple (rank 4 in the current deterministic search):

```bash
python examples/ray/hj_energy_search/audit_whfast_order.py
```

Compare the flagship trajectories with an adaptive IAS15 reference:

```bash
python examples/ray/hj_energy_search/validate_ias15.py
```

Results are written to `best_5/`. `ranking.csv` audits every tested candidate;
each `rank_XX/` directory stores parameters, energy-error samples, and a PNG.
The ranked maximum is evaluated after every fixed step, including steps omitted
from the downsampled CSV and plot.

The hierarchy margin reported in the ranking is

`outer pericentre / largest inner-binary apocentre`.

Only candidates with a margin of at least six are run. This screen is a simple
geometric safeguard, not a mathematical stability proof. The current retained
cases include timestep, particle-order, and IAS15 audits; `best_5/RESULTS.md`
records the results and the longer outer-cycle check.
