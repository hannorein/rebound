if !isdefined(:LIB_REBOUND_INITIALIZED)
const global LIBREBOUND = find_library(["librebound.so"],["/home/ebf11/Code/rebound/","/usr/local/lib","."])
#const rebound_NA_DEFAULT = convert(Cdouble,-2.0); # /* value for transit information that is not determined by TTVFast*/

  # WARNING: rebound has ifdefs that can add additional members and make this type incompatible
  immutable rebound_particle_basic
    x::Cdouble
    y::Cdouble
    z::Cdouble
    vx::Cdouble
    vy::Cdouble
    vz::Cdouble
    ax::Cdouble
    ay::Cdouble
    az::Cdouble
    m::Cdouble
  end

LIB_REBOUND_INITIALIZED = true
end

# librebound

function setp(p::Array{rebound_particle_basic,1})
  ccall( (:setp, LIBREBOUND), Void, (Ptr{rebound_particle_basic},), p )
end

function particle_get(id::Cint)
  retval = ccall( (:particle_get, LIBREBOUND), rebound_particle_basic, (id::Cint,), id )
  return convert(rebound_particle_basic, retval)
  #return particles[i]
end

function particles_get()
  retval = ccall( (:particles_get, LIBREBOUND), Ptr{rebound_particle_basic}, (Void,) )
  return pointer_to_array(retval)
  #	return particles
end

# TODO: Need to relearn how to pass pointers to functions between C and julia
#void set_additional_forces(void (* _cb)(void)){
#	problem_additional_forces = _cb;
#}
#void set_additional_forces_with_parameters(void (* _cb)(struct particle* particles, double t, double dt, double G, int N, int N_megno)){
#	problem_additional_forces_with_parameters = _cb;
#}
#endif

function rebound_step()
  ccall( (:rebound_step, LIBREBOUND), Void, (Void,) )
end

function integrator_set(integrator_id::Cint)
  ccall( (:integrator_set, LIBREBOUND), Void, (Cint,), integrator_id )
end

function rebound_reset()
  ccall( (:reset, LIBREBOUND), Void, () )
end

function integrate!(tmax::Cdouble, exactFinishTime::Cint, keepSynchronized::Cint, particlesModified::Cint, maxR::Cdouble, minD::Cdouble )
  retval = ccall( (:integrate, LIBREBOUND), Cint, (Cdouble, Cint, Cint, Cint, Cdouble, Cdouble),
                 tmax, exactFinishTime, keepSynchronized, particlesModified, maxR, minD)
  return convert(Cint,retval)
end

# particle
get_N() = ccall( (:get_N, LIBREBOUND), CInt, () )
get_Nmax() = ccall( (:get_Nmax, LIBREBOUND), CInt, () )
get_N_active() = ccall( (:get_N_active, LIBREBOUND), CInt, () )
get_N_mengo() = ccall( (:get_N_mengo, LIBREBOUND), CInt, () )
set_Nmax(n::Integer) = ccall( (:set_Nmax, LIBREBOUND), Void, (CInt,), convert(CInt, n) )
set_N_active(n::Integer) = ccall( (:set_N_active, LIBREBOUND), Void, (CInt,), convert(CInt, n) )
set_N_mengo(n::Integer) = ccall( (:set_N_mengo, LIBREBOUND), Void, (CInt,), convert(CInt, n) )
function get_particles()
 p = ccall( (:get_particles, LIBREBOUND), Ptr{Any}, () )
 particles = pointer_to_array(p,1)
 return particles
end


struct particle* particles() { return particles; }

articles_add(p::rebound_particle_basic)  = ccall( (:particles_add, LIBREBOUND), Void, (rebound_particle_basic,), p )
particles_remove_all() = ccall( (:particles_remove_all, LIBREBOUND), Void, () )

# tools
tools_energy() = ccall( (:tools_energy, LIBREBOUND), Cdouble, () )
tools_move_to_center_of_momentum() = ccall( (:tools_move_to_center_of_momentum, LIBREBOUND), Void, () )
tools_init_orbit2d(M::Real, m::Real, a::Real, e::Real, omega::Real, f::Real) =
 ccall( (:tools_init_orbit2d, LIBREBOUND), rebound_particle_basic, (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble), convert(Cdouble,M), convert(Cdouble,m), convert(Cdouble,a), convert(Cdouble,e), convert(Cdouble,omega), convert(Cdouble,f) )
 
tools_init_orbit3d(M::Real, m::Real, a::Real, e::Real, i::Real, Omega::Real, omega::Real, f::Real) =
 ccall( (:tools_init_orbit2d, LIBREBOUND), rebound_particle_basic, (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble), convert(Cdouble,M), convert(Cdouble,m), convert(Cdouble,a), convert(Cdouble,e), convert(Cdouble,i),convert(Cdouble,Omega),convert(Cdouble,omega), convert(Cdouble,f) )

# integrator
integrator_reset() = ccall( (:integrator_reset, LIBREBOUND), Void, () )
integrator_update_acceleration() = ccall( (:integrator_update_acceleration, LIBREBOUND), Void, () )

# gravity_direct
gravity_calculate_acceleration() = ccall( (:gravity_calculate_acceleration, LIBREBOUND), Void, () )
gravity_calculate_variational_acceleration() = ccall( (:gravity_calculate_variational_acceleration, LIBREBOUND), Void, () )

# input
input_binary(filename::Cstring) =  ccall( (:input_binary, LIBREBOUND), Void, (filename,) )

# output  # there are more, but this should get us started
output_ascii(filename::Cstring) =  ccall( (:output_ascii, LIBREBOUND), Void, (filename,) )
output_append_ascii(filename::Cstring) =  ccall( (:output_append_ascii, LIBREBOUND), Void, (filename,) )
output_orbits(filename::Cstring) =  ccall( (:output_orbits, LIBREBOUND), Void, (filename,) )
output_append_orbits(filename::Cstring) =  ccall( (:output_append_orbits, LIBREBOUND), Void, (filename,) )
output_binary(filename::Cstring) =  ccall( (:output_binary, LIBREBOUND), Void, (filename,) )
output_logfile(filename::Cstring) =  ccall( (:output_logfile, LIBREBOUND), Void, (filename,) )
output_png(dirname::Cstring) =  ccall( (:output_png, LIBREBOUND), Void, (dirname,) )
output_png_single(filename::Cstring) =  ccall( (:output_png_single, LIBREBOUND), Void, (filename,) )



function test_rebound_julia()
  rebound_reset()
  # integrator_set(1)

       # dt              = 1e-3; // in year/(2*pi)
       # boxsize         = 3;    // in AU
       # N_active        = 1;    // Only star has non-zero mass. If all particles have mass, delete this line.
        init_box();

        # Initial conditions
        rebound_particle_basic p # // Star
        # The WH integrator assumes a heliocentric coordinate system.
        # Therefore the star has to be at the origin.
        p.x  = 0.0; p.y  = 0.0; p.z  = 0.0;
        p.vx = 0.0; p.vy = 0.0; p.vz = 0.0;
        p.ax = 0.0; p.ay = 0.0; p.az = 0.0;
        p.m  = 1.0; #              // in Solar Masses
        particles_add(p);

        N = 100 # Number of test particles

        eccentricity = 0.4;
        for n=1:N
          p = tools_init_orbit2d(1.0, 0.0, 1.0, eccentricity, 0.,2pi*rand() )
          particles_add(p) # Test particle
        end
end

