# -*- coding: utf-8 -*-
"""Automatically generate citations for a simulation."""

def cite(sim):
    txt = """Simulations in this paper made use of the REBOUND N-body code \citep{rebound}. """
    bib = """@ARTICLE{rebound,
       author = {{Rein}, H. and {Liu}, S. -F.},
        title = "{REBOUND: an open-source multi-purpose N-body code for collisional dynamics}",
      journal = {\\aap},
     keywords = {methods: numerical, planets and satellites: rings, protoplanetary disks, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics, Mathematics - Dynamical Systems, Physics - Computational Physics},
         year = 2012,
        month = jan,
       volume = {537},
          eid = {A128},
        pages = {A128},
          doi = {10.1051/0004-6361/201118085},
archivePrefix = {arXiv},
       eprint = {1110.4876},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2012A&A...537A.128R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""
    if sim.extras:
        txt +="""The REBOUNDx package was used to incorporate additional physics \citep{reboundx}. """
        bib +="""@ARTICLE{reboundx,
       author = {{Tamayo}, Daniel and {Rein}, Hanno and {Shi}, Pengshuai and {Hernandez}, David M.},
        title = "{REBOUNDx: a library for adding conservative and dissipative forces to otherwise symplectic N-body integrations}",
      journal = {\mnras},
     keywords = {gravitation, methods: numerical, planets and satellites: dynamical evolution and stability, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics},
         year = 2020,
        month = jan,
       volume = {491},
       number = {2},
        pages = {2885-2901},
          doi = {10.1093/mnras/stz2870},
archivePrefix = {arXiv},
       eprint = {1908.05634},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""
    if sim.integrator == "ias15":
        txt += """The simulations were integrated using IAS15, a 15th order Gauss-Radau integrator \citep{reboundias15}. """
        bib += """@ARTICLE{reboundias15,
       author = {{Rein}, Hanno and {Spiegel}, David S.},
        title = "{IAS15: a fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits}",
      journal = {\\mnras},
     keywords = {gravitation, methods: numerical, planets and satellites: dynamical evolution and stability, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Solar and Stellar Astrophysics, Mathematics - Numerical Analysis},
         year = 2015,
        month = jan,
       volume = {446},
       number = {2},
        pages = {1424-1437},
          doi = {10.1093/mnras/stu2164},
archivePrefix = {arXiv},
       eprint = {1409.4779},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2015MNRAS.446.1424R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""
    if sim.integrator == "whfast":
        txt += """The simulations were integrated using WHFast, a symplectic Wisdom-Holman integrator \citep{reboundwhfast,wh}. """
        bib += """@ARTICLE{reboundwhfast,
       author = {{Rein}, Hanno and {Tamayo}, Daniel},
        title = "{WHFAST: a fast and unbiased implementation of a symplectic Wisdom-Holman integrator for long-term gravitational simulations}",
      journal = {\\mnras},
     keywords = {gravitation, methods: numerical, planets and satellites: dynamical evolution and stability, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics, Mathematics - Numerical Analysis, Nonlinear Sciences - Chaotic Dynamics, Physics - Computational Physics},
         year = 2015,
        month = sep,
       volume = {452},
       number = {1},
        pages = {376-388},
          doi = {10.1093/mnras/stv1257},
archivePrefix = {arXiv},
       eprint = {1506.01084},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2015MNRAS.452..376R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""
        bib += """@ARTICLE{wh,
       author = {{Wisdom}, Jack and {Holman}, Matthew},
        title = "{Symplectic maps for the N-body problem.}",
      journal = {\aj},
     keywords = {Many Body Problem, Planetary Evolution, Pluto (Planet), Astronomical Maps, Gravitational Effects, Physics (General)},
         year = 1991,
        month = oct,
       volume = {102},
        pages = {1528-1538},
          doi = {10.1086/115978},
       adsurl = {https://ui.adsabs.harvard.edu/abs/1991AJ....102.1528W},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""
        if sim.ri_whfast.kernel != "default":
            txt += """A high order kernel was used in WHFast to improved the accuracy of the integrations \citep{reboundhighorder}. """
            bib += """@ARTICLE{reboundhighorder,
           author = {{Rein}, Hanno and {Tamayo}, Daniel and {Brown}, Garett},
            title = "{High-order symplectic integrators for planetary dynamics and their implementation in REBOUND}",
          journal = {\mnras},
         keywords = {gravitation, methods: numerical, planets and satellites: dynamical evolution and stability, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics, Physics - Computational Physics},
             year = 2019,
            month = nov,
           volume = {489},
           number = {4},
            pages = {4632-4640},
              doi = {10.1093/mnras/stz2503},
    archivePrefix = {arXiv},
           eprint = {1907.11335},
     primaryClass = {astro-ph.EP},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.4632R},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
    """


    
    if sim.integrator == "mercurius":
        txt += """The simulations were integrated using the hybrid symplectic MERCURIUS integrator \citep{reboundmercurius}. """
        bib += """@ARTICLE{reboundmercurius,
       author = {{Rein}, Hanno and {Hernandez}, David M. and {Tamayo}, Daniel and
         {Brown}, Garett and {Eckels}, Emily and {Holmes}, Emma and
         {Lau}, Michelle and {Leblanc}, R{\'e}jean and {Silburt}, Ari},
        title = "{Hybrid symplectic integrators for planetary dynamics}",
      journal = {\mnras},
     keywords = {gravitation, methods: numerical, planets and satellites: dynamical evolution and stability, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics, Mathematics - Dynamical Systems},
         year = 2019,
        month = jun,
       volume = {485},
       number = {4},
        pages = {5490-5497},
          doi = {10.1093/mnras/stz769},
archivePrefix = {arXiv},
       eprint = {1903.04972},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.5490R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""

    if sim.integrator == "janus":
        txt += """The simulations were integrated using the time-reversible JANUS integrator \citep{reboundjanus}. """
        bib += """@ARTICLE{reboundjanus,
       author = {{Rein}, Hanno and {Tamayo}, Daniel},
        title = "{JANUS: a bit-wise reversible integrator for N-body dynamics}",
      journal = {\mnras},
     keywords = {gravitation, methods: numerical, planets and satellites: dynamical evolution and stability, Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Cosmology and Nongalactic Astrophysics, Astrophysics - Earth and Planetary Astrophysics},
         year = 2018,
        month = jan,
       volume = {473},
       number = {3},
        pages = {3351-3357},
          doi = {10.1093/mnras/stx2479},
archivePrefix = {arXiv},
       eprint = {1704.07715},
 primaryClass = {astro-ph.IM},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2018MNRAS.473.3351R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""

    if sim.integrator == "sei":
        txt += """The simulations were integrated using the Symplectic Epicycle Integrator (SEI) \citep{reboundsei}. """
        bib += """@ARTICLE{reboundsei,
       author = {{Rein}, Hanno and {Tremaine}, Scott},
        title = "{Symplectic integrators in the shearing sheet}",
      journal = {\mnras},
     keywords = {methods: numerical, celestial mechanics, planets and satellites: dynamical evolution and stability, planets and satellites: formation, planets and satellites: rings, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Galaxy Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics, Mathematics - Numerical Analysis},
         year = 2011,
        month = aug,
       volume = {415},
       number = {4},
        pages = {3168-3176},
          doi = {10.1111/j.1365-2966.2011.18939.x},
archivePrefix = {arXiv},
       eprint = {1103.1376},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2011MNRAS.415.3168R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""

    if sim.integrator == "saba":
        txt += """The simulations were integrated using the SABA Integrator \citep{reboundhighorder,saba}. """
        bib += """@ARTICLE{reboundhighorder,
       author = {{Rein}, Hanno and {Tamayo}, Daniel and {Brown}, Garett},
        title = "{High-order symplectic integrators for planetary dynamics and their implementation in REBOUND}",
      journal = {\mnras},
     keywords = {gravitation, methods: numerical, planets and satellites: dynamical evolution and stability, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics, Physics - Computational Physics},
         year = 2019,
        month = nov,
       volume = {489},
       number = {4},
        pages = {4632-4640},
          doi = {10.1093/mnras/stz2503},
archivePrefix = {arXiv},
       eprint = {1907.11335},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.4632R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""
        bib += """@ARTICLE{saba,
       author = {{Laskar}, Jacques and {Robutel}, Philippe},
        title = "{High order symplectic integrators for perturbed Hamiltonian systems}",
      journal = {Celestial Mechanics and Dynamical Astronomy},
     keywords = {SYMPLECTIC INTEGRATORS, HAMILTONIAN SYSTEMS, PLANETARY MOTION, LIE ALGEBRA, Astrophysics},
         year = 2001,
        month = jul,
       volume = {80},
       number = {1},
        pages = {39-62},
archivePrefix = {arXiv},
       eprint = {astro-ph/0005074},
 primaryClass = {astro-ph},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2001CeMDA..80...39L},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""


    if sim.simulationarchive_auto_interval!=0 or sim.simulationarchive_auto_walltime!=0 or sim.simulationarchive_auto_step!=0:
        txt += """The SimulationArchive format was used to store fully reproducible simulation data \citep{reboundsa}. """
        bib += """@ARTICLE{reboundsa,
       author = {{Rein}, Hanno and {Tamayo}, Daniel},
        title = "{A new paradigm for reproducing and analyzing N-body simulations of planetary systems}",
      journal = {\mnras},
     keywords = {methods: numerical, gravitation, planets and satellites: dynamical evolution and stability, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics},
         year = 2017,
        month = may,
       volume = {467},
       number = {2},
        pages = {2377-2383},
          doi = {10.1093/mnras/stx232},
archivePrefix = {arXiv},
       eprint = {1701.07423},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.2377R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""
    if sim.N_var>0:
        txt += """Variational equations were used to calculate trajectories of nearby orbits \citep{reboundvar}. """
        bib += """@ARTICLE{reboundvar,
       author = {{Rein}, Hanno and {Tamayo}, Daniel},
        title = "{Second-order variational equations for N-body simulations}",
      journal = {\mnras},
     keywords = {gravitation, methods: numerical, planets and satellites: dynamical evolution and stability, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics, Mathematics - Classical Analysis and ODEs, Mathematics - Dynamical Systems},
         year = 2016,
        month = jul,
       volume = {459},
       number = {3},
        pages = {2275-2285},
          doi = {10.1093/mnras/stw644},
archivePrefix = {arXiv},
       eprint = {1603.03424},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2016MNRAS.459.2275R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""




    return txt, bib
