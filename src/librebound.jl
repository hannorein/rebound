if !isdefined(:LIB_REBOUND_INITIALIZED)
const global LIBREBOUND = find_library(["librebound.so"],["/home/ebf11/Code/rebound/","/usr/local/lib",".","../../"])

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

function rebound_particle_basic(_x::Real,_y::Real,_z::Real, _vx::Real,_vy::Real,_vz::Real; _m::Real = 0.0)
  return rebound_particle_basic(_x,_y,_z,_vx,_vy,_vz,0.0,0.0,0.0,_m)
end

function rebound_particle_basic(_m::Real)
  return rebound_particle_basic(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,_m)
end

# librebound

function setp(p::Array{rebound_particle_basic,1})
  ccall( (:setp, LIBREBOUND), Void, (Ptr{rebound_particle_basic},), p )
end

# Account for different indexing
function particle_get(id::Integer)
  @assert(1<=id<=get_N())
  retval = ccall( (:particle_get, LIBREBOUND), rebound_particle_basic, (Cint,), convert(Cint,id-1) )
  return convert(rebound_particle_basic, retval)
end

function particles_get()
  retval = ccall( (:particles_get, LIBREBOUND), Ptr{rebound_particle_basic}, (Void,) )
  return pointer_to_array(retval)
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

# Should do something so Julia doesn't need to use magic numbers here
function integrator_set(integrator_id::Integer)
  @assert(0<=integrator_id<=6)
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

# particle
# accesor/mutator functions
get_N() = ccall( (:get_N, LIBREBOUND), Cint, () )
get_Nmax() = ccall( (:get_Nmax, LIBREBOUND), Cint, () )
get_N_active() = ccall( (:get_N_active, LIBREBOUND), Cint, () )
get_N_mengo() = ccall( (:get_N_mengo, LIBREBOUND), Cint, () )
set_N(n::Integer) = ccall( (:set_N, LIBREBOUND), Void, (Cint,), convert(Cint, n) )
set_Nmax(n::Integer) = ccall( (:set_Nmax, LIBREBOUND), Void, (Cint,), convert(Cint, n) )
set_N_active(n::Integer) = ccall( (:set_N_active, LIBREBOUND), Void, (Cint,), convert(Cint, n) )
set_N_mengo(n::Integer) = ccall( (:set_N_mengo, LIBREBOUND), Void, (Cint,), convert(Cint, n) )

function get_particles()
 num_part = get_N()
 #p = ccall( (:get_particles, LIBREBOUND), Ptr{Any}, () )
 #particles = pointer_to_array(convert(Ptr{rebound_particle_basic},p),num_part)
 p = ccall( (:get_particles, LIBREBOUND), Ptr{rebound_particle_basic}, () )
 particles = pointer_to_array(p,num_part)
 return particles
end

# WARNING: Trying direct/unsafe access w/o access functions, before deciding whether to use these instead of C accessor functions
get_N_direct() = unsafe_load(convert(Ptr{Cint},cglobal((:N,LIBREBOUND))))
set_N_direct(n::Integer) = unsafe_store!(convert(Ptr{Cint},cglobal((:N,LIBREBOUND))),convert(Cint,n))

get_dt() = unsafe_load(convert(Ptr{Cdouble},cglobal((:dt,LIBREBOUND))))
set_dt(dt::Real) = unsafe_store!(convert(Ptr{Cdouble},cglobal((:dt,LIBREBOUND))),convert(Cdouble,dt))

particles_add(p::rebound_particle_basic)  = ccall( (:particles_add_jl, LIBREBOUND), Void, (Ptr{rebound_particle_basic},), &p )
particles_remove_all() = ccall( (:particles_remove_all, LIBREBOUND), Void, () )

# WARNING:  Most everything below here is totally untested
# tools
tools_energy() = ccall( (:tools_energy, LIBREBOUND), Cdouble, () )
tools_move_to_center_of_momentum() = ccall( (:tools_move_to_center_of_momentum, LIBREBOUND), Void, () )
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
output_png(dirname::Cstring) =  ccall( (:output_png, LIBREBOUND), Void, (dirname,) )
output_png_single(filename::Cstring) =  ccall( (:output_png_single, LIBREBOUND), Void, (filename,) )


