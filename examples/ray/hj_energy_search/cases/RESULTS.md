# Current fixed-hierarchy search results

## Result

A stable 2+2 quadruple provides a clear case where the fixed hierarchy

`[[A1,A2],[B1,B2]]`

is dramatically better than ordinary WHFast's linear Jacobi comb. The search
evaluated 640 deterministic parameter combinations. The score is

`max relative energy error of WHFast / max relative energy error of fixed HJ`.

The maximum is evaluated after every fixed step, not only at plot samples.

| Rank | Candidate | WHFast max error | Fixed-HJ max error | Search ratio |
|---:|---|---:|---:|---:|
| 1 | `candidate_0169` | 5.562e-3 | 8.460e-10 | 6.574e6 |
| 2 | `candidate_0423` | 1.839e-2 | 7.253e-9 | 2.535e6 |
| 3 | `candidate_0210` | 4.857e-3 | 2.388e-9 | 2.034e6 |
| 4 | `candidate_0155` | 1.464e-2 | 8.174e-9 | 1.791e6 |
| 5 | `candidate_0375` | 1.739e-3 | 1.513e-9 | 1.149e6 |

## Flagship configuration

`candidate_0169` consists of two tight binaries on a wide mutual orbit:

- masses `(1.3, 0.3)` and `(0.7, 0.7)`;
- inner semimajor axes `(1.0, 0.55)` and eccentricities `(0.5, 0.4)`;
- outer semimajor axis `80`, eccentricity `0.5`;
- mutual inclination `145 degrees`;
- initial true anomalies `(73, 211, 179) degrees`;
- timestep `P_ref/40`;
- outer-pericentre / largest-inner-apocentre margin `26.67`.

The original 120-slow-inner-period run gives a 6.57-million-fold energy-error
ratio. A five-times-longer run spans about 1.15 outer periods, including the
wide orbit's pericentre region. On that longer run the maximum errors are
`5.568e-3` for WHFast and `2.670e-8` for fixed HJ, a `2.086e5` ratio.

## Robustness checks

- Timestep scan: at `P_ref/{20,40,80,160}`, the flagship's longer-run ratio is
  between `1.87e5` and `2.11e5`; both methods show bounded error and the errors
  decrease as the timestep decreases. The deliberately coarse `P_ref/10` run
  is outside the 10,000-fold target on the full outer-cycle test.
- Particle ordering: all 24 ordinary-WHFast orders were tested. Even its best
  order has error `7.596e-4`, still `8.98e5` times the fixed-HJ search error.
- IAS15: with `epsilon=1e-12` over the five-times-longer run, IAS15's sampled
  maximum energy error is `8.28e-14`. Relative to IAS15, ordinary WHFast's
  hierarchy-coordinate RMS errors are `8.93e-3` in position and `1.11` in
  velocity; fixed HJ's are `3.98e-10` and `5.94e-8`, respectively.
- Every current top-five case exceeds 10,000-fold at its retained timestep,
  across the longer run, and against the best of all 24 WHFast orderings.

## Why this works

Ordinary Jacobi coordinates form a chain. They can make one tight binary an
exact Kepler subproblem, but cannot simultaneously encode two disjoint tight
binaries. The branching HJ tree makes both inner binaries and their wide mutual
orbit explicit Kepler nodes, leaving only the genuinely small non-Keplerian
interaction in the kick. The advantage is therefore structural rather than a
close-encounter effect.

This agrees with Beust's hierarchical-Jacobi formulation for multiple stellar
systems and with later work using hierarchical Jacobi coordinates for triple
stars:

- Beust (2003), *Symplectic integration of hierarchical stellar systems*:
  https://www.aanda.org/articles/aa/pdf/2003/12/aa3133.pdf
- Verrier & Evans (2007), *Planetary Stability Zones in Hierarchical Triple
  Star Systems*: https://arxiv.org/abs/0710.1167
- REBOUND WHFast coordinate documentation:
  https://rebound.hanno-rein.de/ipython_examples/AdvWHFast/

These are constructed numerical benchmarks, not observed systems. Their value
is to isolate the Hamiltonian-splitting advantage of the correct fixed tree.
