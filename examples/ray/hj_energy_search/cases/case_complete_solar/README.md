# case_complete_solar

This benchmark asks whether a fixed hierarchical-Jacobi tree improves a
Solar-System model when several major moons are resolved as massive particles.

It contains the Sun; Mercury through Neptune; Pluto; the Moon; Io, Europa,
Ganymede, and Callisto; Titan; Titania; Triton; and Charon. Approximate masses
and orbital scales are used with deterministic phases. It is a controlled
dynamical benchmark rather than an epoch-matched JPL ephemeris.

The physical fixed tree is a heliocentric radial comb containing these local
branches:

```text
[Earth,Moon]
[[[[Jupiter,Io],Europa],Ganymede],Callisto]
[Saturn,Titan]
[Uranus,Titania]
[Neptune,Triton]
[Pluto,Charon]
```

Run the default five-year comparison and timestep scan:

```bash
python examples/ray/hj_energy_search/cases/case_complete_solar/run.py
```

The script compares ordinary WHFast and fixed-tree WHFast-HJ at identical
Cartesian initial conditions and timesteps, checks trajectory errors against
IAS15, and writes CSV, JSON, and PNG results in this directory.

## Current five-year result

The timestep is normalized to Io's 1.76925-day orbital period. At the main
choice `dt = P_Io/20`:

```text
WHFast maximum total relative-energy error:  5.711e-8
fixed-HJ maximum total relative-energy error: 4.284e-11
WHFast / HJ:                                 1.333e3
IAS15 sampled maximum error:                 7.042e-15
```

The IAS15-referenced full-system RMS trajectory errors are:

```text
                         position       velocity
WHFast                   1.509e-4       3.944e-1
fixed HJ                 6.459e-9       1.064e-5
```

HJ is therefore clearly better for this resolved-moon model, although the
total-energy improvement is about 1,300 times rather than the million-fold
improvement seen in the synthetic equal-strength 2+2 quadruples. At very small
timesteps the HJ energy curve approaches a roundoff/accumulation floor, so the
ratio falls even though HJ remains more accurate.

This case is deliberately more realistic than the synthetic search cases, but
it is not a claim about a production Solar-System ephemeris: the orbital phases
are constructed, relativistic precession and oblateness are omitted, and only
selected major moons are resolved.
