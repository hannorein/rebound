if !isdefined(:LIB_REBOUND_INITIALIZED)
  const global LIBREBOUND = find_library(["librebound.so"],["/home/ebf11/Code/rebound/","/usr/local/lib",".","../../"])
  
  # REBOUND particle structure
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

  # Integrator lookup dictionary
  INTEGRATORS = {"ias15" => 0, "whfast" => 1, "sei" => 2, "wh" => 3, "leapfrog" => 4, "hybrid" => 5, "none" => 6}

  LIB_REBOUND_INITIALIZED = true
end

function rebound_particle_basic(_x::Real,_y::Real,_z::Real, _vx::Real,_vy::Real,_vz::Real; _m::Real = 0.0)
  return rebound_particle_basic(_x,_y,_z,_vx,_vy,_vz,0.0,0.0,0.0,_m)
end

function rebound_particle_basic(_m::Real)
  return rebound_particle_basic(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,_m)
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

function set_integrator(integrator_name::ASCIIString)
  local integrator_id = INTEGRATORS[integrator_name]
  ccall( (:integrator_set, LIBREBOUND), Void, (Cint,), convert(Cint,integrator_id) )
end

function rebound_reset()
  ccall( (:reset, LIBREBOUND), Void, () )
end

function integrate(tmax::Real, exactFinishTime::Integer=0, keepSynchronized::Integer=0, particlesModified::Integer=0, maxR::Real=0.0, minD::Real=0.0 )
  retval = ccall( (:integrate, LIBREBOUND), Cint, 
		(Cdouble, Cint, Cint, Cint, Cdouble, Cdouble),
                 convert(Cdouble,tmax), convert(Cint,exactFinishTime), convert(Cint,keepSynchronized), convert(Cint,particlesModified), convert(Cdouble,maxR), convert(Cdouble,minD) )
  return convert(Cint,retval)
end

get_N() = unsafe_load(convert(Ptr{Cint},cglobal((:N,LIBREBOUND))))
set_N(n::Integer) = unsafe_store!(convert(Ptr{Cint},cglobal((:N,LIBREBOUND))),convert(Cint,n))
set_N_megnopp(n::Integer) = unsafe_store!(convert(Ptr{Cint},cglobal((:N_megnopp,LIBREBOUND))),convert(Cint,n))

get_dt() = unsafe_load(convert(Ptr{Cdouble},cglobal((:dt,LIBREBOUND))))
set_dt(dt::Real) = unsafe_store!(convert(Ptr{Cdouble},cglobal((:dt,LIBREBOUND))),convert(Cdouble,dt))

function get_particles()
 num_part = get_N()
 p = unsafe_load(cglobal((:particles,LIBREBOUND), Ptr{rebound_particle_basic}))
 particles = pointer_to_array(p,num_part)
 return particles
end

function get_particle(id::Integer)
  return get_particles(id)
end

add(p::rebound_particle_basic)  = ccall( (:particles_add_ptr, LIBREBOUND), Void, (Ptr{rebound_particle_basic},), &p )
particles_remove_all() = ccall( (:particles_remove_all, LIBREBOUND), Void, () )

# WARNING:  Most everything below here is totally untested
# tools
tools_energy() = ccall( (:tools_energy, LIBREBOUND), Cdouble, () )
tools_move_to_com() = ccall( (:tools_move_to_center_of_momentum, LIBREBOUND), Void, () )
tools_init_orbit2d(M::Real, m::Real, a::Real, e::Real, omega::Real, f::Real) =
 ccall( (:tools_init_orbit2d, LIBREBOUND), rebound_particle_basic, (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble), convert(Cdouble,M), convert(Cdouble,m), convert(Cdouble,a), convert(Cdouble,e), convert(Cdouble,omega), convert(Cdouble,f) )
 
tools_init_orbit3d(M::Real, m::Real, a::Real, e::Real, i::Real, Omega::Real, omega::Real, f::Real) =
 ccall( (:tools_init_orbit2d, LIBREBOUND), rebound_particle_basic, (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble), convert(Cdouble,M), convert(Cdouble,m), convert(Cdouble,a), convert(Cdouble,e), convert(Cdouble,i),convert(Cdouble,Omega),convert(Cdouble,omega), convert(Cdouble,f) )

tools_megno_init(delta::Real) = ccall( (:tool_megno_init, LIBREBOUND), Void, (Cdouble,), convert(Cdoubl, delta) )


# integrator
integrator_reset() = ccall( (:integrator_reset, LIBREBOUND), Void, () )
integrator_update_acceleration() = ccall( (:integrator_update_acceleration, LIBREBOUND), Void, () )

# gravity_direct
gravity_calculate_acceleration() = ccall( (:gravity_calculate_acceleration, LIBREBOUND), Void, () )
gravity_calculate_variational_acceleration() = ccall( (:gravity_calculate_variational_acceleration, LIBREBOUND), Void, () )

# input
# WARNING: Figure out how to check julia version, for compatability with v0.4
typealias Cstring Ptr{Uint8}

input_binary(filename::Cstring) =  ccall( (:input_binary, LIBREBOUND), Void, (filename,) )

# output  # there are more, but this should get us started
output_ascii(filename::Cstring) =  ccall( (:output_ascii, LIBREBOUND), Void, (filename,) )
output_append_ascii(filename::Cstring) =  ccall( (:output_append_ascii, LIBREBOUND), Void, (filename,) )
output_orbits(filename::Cstring) =  ccall( (:output_orbits, LIBREBOUND), Void, (filename,) )
output_append_orbits(filename::Cstring) =  ccall( (:output_append_orbits, LIBREBOUND), Void, (filename,) )
output_binary(filename::Cstring) =  ccall( (:output_binary, LIBREBOUND), Void, (filename,) )
output_logfile(filename::Cstring) =  ccall( (:output_logfile, LIBREBOUND), Void, (filename,) )


