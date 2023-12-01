import ctypes
from . import GenericError
from .particle import Particle
from .simulation import Simulation

class Variation(ctypes.Structure):
    """
    REBOUND Variational Configuration Object.

    This object encapsulated the configuration of one set of variational 
    equations in a REBOUND simulation.  It is an abstraction of the 
    C struct reb_variational_configuration.

    None of the fields in this struct should be changed after it has
    been initialized.

    One rebound simulation can include any number of first and second order 
    variational equations.

    Note that variations are only encoded as particles for convenience.  
    A variational particle's position and velocity should be interpreted as a derivative, i.e. how much that position or velocity varies with respect to the first or second-order variation.  
    See ipython_examples/VariationalEquations.ipynb and Rein and Tamayo (2016) for details.

    Examples
    --------

    >>> sim = rebound.Simulation()          # Create a simulation
    >>> sim.add(m=1.)                       # Add a star
    >>> sim.add(m=1.e-3, a=1.)              #     a planet
    >>> var_config = sim.add_variation()    # Add a set of variational equations. 
    >>> var_config.particles[1].x = 1.      # Initialize the variational particle corresponding to the planet
    >>> sim.integrate(100.)                 # Integrate the simulation
    >>> print(var_config.particles[0].vx)   # Print the velocity of the variational particle corresponding to the star
    """
    _fields_ = [
                ("_sim", ctypes.POINTER(Simulation)),
                ("order", ctypes.c_int),
                ("index", ctypes.c_int),
                ("testparticle", ctypes.c_int),
                ("index_1st_order_a", ctypes.c_int),
                ("index_1st_order_b", ctypes.c_int),
                ("_lrescale", ctypes.c_double)]

    def vary(self, particle_index, variation, variation2=None, primary=None):
        """
        This function can be used to initialize the variational particles that are 
        part of a Variation.
    
        Note that rather than using this convenience function, one can 
        also directly manipulate the particles' coordinate using the following
        syntax:

        >>> var = sim.add_variation()
        >>> var.particles[0].x = 1.

        The ``vary()`` function is useful for initializing variations corresponding to 
        changes in one of the orbital parameters for a particle on a bound 
        Keplerian orbit.

        The function supports both first and second order variations in the following
        classical orbital parameters:
          a, e, inc, omega, Omega, f
        as well as the Pal (2009) coordinates: 
          a, h, k, ix, iy, lambda
        and in both cases the mass m of the particle. The advantage of the Pal coordinate
        system is that all derivatives are well behaved (infinitely differentiable).
        Classical orbital parameters on the other hand exhibit coordinate singularities, 
        for example when e=0.
        
        The following example initializes the variational particles corresponding to a 
        change in the semi-major axis of the particle with index 1:
        
        >>> var = sim.add_variation()
        >>> var.vary(1,"a")

        Parameters
        ----------
        particle_index : int
            The index of the particle that should be varied. The index starts at 0 and runs through N-1. The first particle added to a simulation receives the index 0, the second 1, and the on.
        variation : string
            This parameter determines which orbital parameter is varied. 
        variation2: string, optional
            This is only used for second order variations which can depend on two varying parameters. If omitted, then it is assumed that the parameter variation is variation2.
        primary: Particle, optional
            By default, variational particles are created in the Heliocentric frame. 
            Set this parameter to use any other particles as a primary (e.g. the center of mass).
        """
        if self.order==2 and variation2 is None:
            variation2 = variation
        if self._sim is not None:
            sim = self._sim.contents
            particles = sim.particles
        else:
            raise RuntimeError("Something went wrong. Cannot seem to find simulation corresponding to variation.")
        if self.testparticle >= 0:
            particles[self.index] = Particle(simulation=sim,particle=particles[particle_index], variation=variation, variation2=variation2, primary=primary)
        else:
            particles[self.index + particle_index] = Particle(simulation=sim,particle=particles[particle_index], variation=variation, variation2=variation2, primary=primary)

    @property
    def particles(self):
        """
        Access the variational particles corresponding to this set of variational equations.

        The function returns a list of particles which are sorted in the same way as those in 
        sim.particles

        The particles are pointers and thus can be modified. 

        If there are N real particles, this function will also return a list of N particles (all of which are 
        variational particles). 
        """
        sim = self._sim.contents
        ps = []
        if self.testparticle>=0:
            N = 1
        else:
            N = sim.N-sim.N_var 
        
        ParticleList = Particle*N
        ps = ParticleList.from_address(ctypes.addressof(sim._particles.contents)+self.index*ctypes.sizeof(Particle))
        return ps
    
    @property
    def lrescale(self):
        """
        Access the lrescale parameter. 
        
        This is a property because sim.add_variation() returns a copy of the struct, so need to find up-to-date reb_variational_configuration struct in simulation.
        """
        sim = self._sim.contents
        for i in range(sim.N_var_config):
            if sim.var_config[i].index == self.index:
                return sim.var_config[i]._lrescale
        raise GenericError("An error occured while trying to find variational struct in simulation.")
    
    @lrescale.setter
    def lrescale(self, value):
        """
        Set the lrescale parameter. 
        
        This is a property because sim.add_variation() returns a copy of the struct, so need to find up-to-date reb_variational_configuration struct in simulation.
        """
        sim = self._sim.contents
        for i in range(sim.N_var_config):
            if sim.var_config[i].index == self.index:
                self._lrescale = ctypes.c_double(value)
                sim.var_config[i]._lrescale = ctypes.c_double(value)
                return

        raise GenericError("An error occured while trying to find variational struct in simulation.")

