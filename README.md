[![Version](https://img.shields.io/badge/rebound-v4.4.10-green.svg?style=flat)](https://rebound.readthedocs.org)
[![codecov](https://codecov.io/github/hannorein/rebound/graph/badge.svg?token=Zmynoi99Vl)](https://codecov.io/github/hannorein/rebound)
[![PyPI](https://badge.fury.io/py/rebound.svg)](https://badge.fury.io/py/rebound)
[![GPL](https://img.shields.io/badge/license-GPL-green.svg?style=flat)](https://github.com/hannorein/rebound/blob/main/LICENSE)
[![Paper](https://img.shields.io/badge/arXiv-1110.4876-green.svg?style=flat)](https://arxiv.org/abs/1110.4876)
[![Paper](https://img.shields.io/badge/arXiv-1409.4779-green.svg?style=flat)](https://arxiv.org/abs/1409.4779)
[![Paper](https://img.shields.io/badge/arXiv-1506.01084-green.svg?style=flat)](https://arxiv.org/abs/1506.01084)
[![Paper](https://img.shields.io/badge/arXiv-1603.03424-green.svg?style=flat)](https://arxiv.org/abs/1603.03424)
[![Paper](https://img.shields.io/badge/arXiv-1701.07423-green.svg?style=flat)](https://arxiv.org/abs/1701.07423)
[![Paper](https://img.shields.io/badge/arXiv-1704.07715-green.svg?style=flat)](https://arxiv.org/abs/1704.07715)
[![Paper](https://img.shields.io/badge/arXiv-1903.04972-green.svg?style=flat)](https://arxiv.org/abs/1903.04972)
[![Paper](https://img.shields.io/badge/arXiv-1907.11335-green.svg?style=flat)](https://arxiv.org/abs/1907.11335)
[![Docs](https://readthedocs.org/projects/rebound/badge/?version=latest)](https://rebound.readthedocs.io/en/latest/?badge=latest)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/hannorein/rebound/main)
[![REBOUND (C)](https://github.com/hannorein/rebound/actions/workflows/c.yml/badge.svg)](https://github.com/hannorein/rebound/actions/workflows/c.yml)
[![REBOUND (python)](https://github.com/hannorein/rebound/actions/workflows/python.yml/badge.svg)](https://github.com/hannorein/rebound/actions/workflows/python.yml)
    

# Welcome to REBOUND

![REBOUND Examples](https://github.com/hannorein/rebound/raw/main/docs/img/reboundbanner.png)

REBOUND is an N-body integrator, i.e. a software package that can integrate the motion of particles under the influence of gravity. The particles can represent stars, planets, moons, ring or dust particles. REBOUND is very flexible and can be customized to accurately and efficiently solve many problems in astrophysics.  

## Features

* No dependencies on external libraries.
* Runs natively on Linux, MacOS, and Windows. 
* Symplectic integrators WHFast, SEI, LEAPFROG, EOS.
* Hybrid symplectic integrators for planetary dynamics with close encounters MERCURIUS
* Hybrid reversible integrators for planetary dynamics with arbitrary close encounters TRACE
* High order symplectic integrators for integrating planetary systems SABA, WH Kernel methods.
* High accuracy non-symplectic integrator with adaptive time-stepping IAS15.
* Can integrate arbitrary user-defined ODEs that are coupled to N-body dynamics for tides, spin, etc
* Support for collisional/granular dynamics, various collision detection routines
* The computationally intensive parts of the code are written entirely in C, conforming to the ISO standard C99, and can be used as a thread-safe shared library
* Easy-to-use Python module, installation in 3 words: `pip install rebound`
* Real-time, 3D visualization, for both C and Python.
* Extensive set of example problems for both C and Python. You can run examples directly from your browser without the need to download or install anything.
* Parallelized WHFast512 integrator for super fast integrations of planetary systems with SIMD AVX512 instructions
* Parallelized with OpenMP (for shared memory systems)
* Parallelized with MPI is supported for some special use cases only (using an essential tree for gravity and collisions)
* The code is 100% open-source. All features are included in the public repository on github.

## Try out REBOUND 

You can try out REBOUND without installing it. 
Simply head over to [readthedocs.org](https://rebound.readthedocs.io/en/latest/examples/).
All the C examples have been compiled with emscripten and can run directly in your browser.

## One minute installation

You can install REBOUND with pip if you want to only use the python version of REBOUND:

    pip install rebound

Then, you can run a simple REBOUND simulation such as

```python
import rebound
sim = rebound.Simulation()
sim.add(m=1.0)
sim.add(m=1.0e-3, a=1.0)
sim.integrate(1000.)
sim.status()
```

If you want to use the C version of REBOUND simply copy and paste this line into your terminal (it won't do anything bad, we promise):

```bash
git clone https://github.com/hannorein/rebound && cd rebound/examples/shearing_sheet && make && ./rebound
```

 
## Documentation
The full documentation with many examples, changelogs and tutorials can be found at

<https://rebound.readthedocs.org>

If you have trouble installing or using REBOUND, please open an issue on github and we'll try to help as much as we can.

There are also short YouTube videos describing various aspects of REBOUND available at https://www.youtube.com/channel/UCNmrCzxcmWVTBwtDPPLxkkw .

## Related projects

### Additional physics
To easily incorporate additional physics modules such as migration forces, GR effects and spin into your REBOUND simulations, see REBOUNDx at https://github.com/dtamayo/reboundx

### Analytical and semianalytical tools
If you're interested in comparing numerical simulations to analytical and semianalytical tools for celestial mechanics, see Celmech at https://github.com/shadden/celmech

### Ephemeris-quality integrations of test particles
To generate ephemeris-quality integrations of test particles in the Solar System with a precision on par with JPL's small body integrator, see ASSIST at https://github.com/matthewholman/assist

## Papers

There are several papers describing the functionality of REBOUND.

1. Rein & Liu 2012 (Astronomy and Astrophysics, Volume 537, A128) describes the code structure and the main feature including the gravity and collision routines for many particle systems. <http://adsabs.harvard.edu/abs/2012A%26A...537A.128R>

2. Rein & Tremaine 2011 (Monthly Notices of the Royal Astronomical Society, Volume 415, Issue 4, pp. 3168-3176) describes the Symplectic Epicycle integrator for shearing sheet simulations. <https://ui.adsabs.harvard.edu/abs/2011MNRAS.415.3168R>

3. Rein & Spiegel 2015 (Monthly Notices of the Royal Astronomical Society, Volume 446, Issue 2, p.1424-1437) describes the versatile high order integrator IAS15 which is now part of REBOUND. <http://adsabs.harvard.edu/abs/2015MNRAS.446.1424R>

4. Rein & Tamayo 2015 (Monthly Notices of the Royal Astronomical Society, Volume 452, Issue 1, p.376-388) describes WHFast, the fast and unbiased implementation of a symplectic Wisdom-Holman integrator for long term gravitational simulations. <http://adsabs.harvard.edu/abs/2015MNRAS.452..376R>

5. Rein & Tamayo 2016 (Monthly Notices of the Royal Astronomical Society, Volume 459, Issue 3, p.2275-2285) develop the framework for second order variational equations. <https://ui.adsabs.harvard.edu/abs/2016MNRAS.459.2275R>

6. Rein & Tamayo 2017 (Monthly Notices of the Royal Astronomical Society, Volume 467, Issue 2, p.2377-2383) describes the Simulationarchive for exact reproducibility of N-body simulations. <https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.2377R>

7. Rein & Tamayo 2018 (Monthly Notices of the Royal Astronomical Society, Volume 473, Issue 3, p.3351–3357) describes the integer based JANUS integrator. <https://ui.adsabs.harvard.edu/abs/2018MNRAS.473.3351R>

8. Rein, Hernandez, Tamayo, Brown, Eckels, Holmes, Lau, Leblanc & Silburt 2019 (Monthly Notices of the Royal Astronomical Society, Volume 485, Issue 4, p.5490-5497) describes the hybrid symplectic integrator MERCURIUS. <https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.5490R>

9. Rein, Tamayo & Brown 2019 (Monthly Notices of the Royal Astronomical Society, Volume 489, Issue 4, November 2019, Pages 4632-4640) describes the implementation of the high order symplectic integrators SABA, SABAC, SABACL, WHCKL, WHCKM, and WHCKC. <https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.4632R/>

## Acknowledgments

If you use this code or parts of this code for results presented in a scientific publication, we would greatly appreciate a citation.
The simplest way to find the citations relevant to the specific setup of your REBOUND simulation is: 

```python
sim = rebound.Simulation()
-your setup-
sim.cite()
```


## Contributors

* Hanno Rein, University of Toronto, <hanno@hanno-rein.de>
* Dan Tamayo, Harvey Mudd College, <dtamayo@hmc.edu>
* David S. Spiegel, Institute for Advanced Study Princeton, <dave@ias.edu>
* Garett Brown, University of Toronto, <garett.brown@mail.utoronto.ca>
* Shangfei Liu, Kavli Institute for Astronomy and Astrophysics at Peking University, <liushangfei@pku.edu.cn>
* Ari Silburt, Penn State University, <ajs725@psu.edu>
* and many others! Check the git history to find out who contributed to the code.

REBOUND is open source and you are invited to contribute to this project! 


## License

REBOUND is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

REBOUND is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with REBOUND.  If not, see <http://www.gnu.org/licenses/>.

