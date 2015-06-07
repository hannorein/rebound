include("../../src/librebound.jl")

function test_rebound_julia()
       rebound_reset()
       integrator_set(1)
       set_dt(1.e-3)        # in year/(2*pi)
       # WARNING:  I don't know if we need to set these
       # boxsize         = 3;    // in AU
       # N_active        = 1;    // Only star has non-zero mass. If all particles have mass, delete this line.
       # init_box();

        # Initial conditions
	star_mass = 1.0
        pstar = rebound_particle_basic(star_mass) # // Star
        particles_add(pstar);  # The WH integrator assumes a heliocentric coordinate system.  # Therefore the star has to be at the origin.

        local N = 10 # Number of particles
        for n=1:N
          local m_pl = 1e-4*rand()
          local a = 1.0
          local anom = 2pi*rand()
          local eccentricity = 0.4*rand()
          local omega = 2pi*rand()
          p = tools_init_orbit2d(star_mass, m_pl, a, eccentricity, omega,anom )
          particles_add(p) # Test particle
        end
    local e_init = tools_energy()
    println("# Starting integrateion: E=$e_init")
    local t_int = 10.0*2pi
    local flag = integrate(t_int)
    local e_final = tools_energy()
    println("# Ending integrateion: E=$e_init")
    local delta_e = e_final - e_init
    println("# Integrated for $t_int, flag=$flag, dE/E = ",delta_e/e_init)
end

test_rebound_julia()

