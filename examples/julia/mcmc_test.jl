include("librebound.jl")
using Distributions
using Lora
# using PyPlot

# Simulated Data:
data_N = 30
data_t = linspace(0.1,100,data_N).+0.5.*rand(data_N)
data_rv = 0.500883379885605.*(sin(2pi.*data_t/50.+0.2*pi).+0.6.*(rand(data_N).-0.5))

# plot(data_t,data_rv)

# Setting up simulations
function setup_rv_simulation( param::Array{Float64,1} )
  @assert(length(param)==4)
  star_mass = 1.0
  a=param[1]
  ecc=param[2]
  omega=param[3]
  anom=param[4]
  m_pl = 1e-3
  rebound_reset()
  # integrator_set(1)
  # set_dt(1.e-3)        # in year/(2*pi)

  add( rebound_particle_basic(1.0) )
  add( tools_init_orbit2d(star_mass, m_pl, a, ecc, omega, anom) )
end

# Calc RVs
function calc_rvs( data_t::Array{Float64,1} )
  @assert(length(data_t)>=1)
  ps = get_particles()
  sim_rv = zeros(length(data_t))
  for (i,t) in enumerate(data_t)
     integrate(t)
     sim_rv[i] = ps[2].vx
  end
  return sim_rv
end

# Calc RVs, chisq and its gradient
function calc_rvs_chisq_grad( param::Array{Float64,1} )
  @assert(length(param)==4)
  star_mass = 1.0
  # a=param[1]; ecc=param[2]; omega=param[3]; anom=param[4]
  m_pl = 1e-3
  rebound_reset()
  set_N_megnopp(4)
  tools_megno_init(0.)
  # integrator_set(1)
  # set_dt(1.e-3)        # in year/(2*pi)

  add( rebound_particle_basic(star_mass) )
  add( tools_init_orbit2d(star_mass, m_pl, param[1], param[2], param[3], param[4]) )
  ps = get_particles()
  N = get_N() # this is the total number of particles (incl variational particles)
  for k in [1:4]
    delta = 1e-6
	param2 = param
	param2[k] += delta

	mp = tools_init_orbit2d(star_mass, m_pl, param2[1], param2[2], param2[3], param2[4]) 
	vari = (N*k)+2
	ps[vari].x = ps[1].x - mp.x
	ps[vari].y = ps[1].y - mp.y
	ps[vari].z = ps[1].z - mp.z
	ps[vari].vx = ps[1].vx - mp.vx
	ps[vari].vy = ps[1].vy - mp.vy
	ps[vari].vz = ps[1].vz - mp.vz
  end

  sim_rv = zeros(length(data_t))
  chisq = 0.0
  d_chisq = 0.0
  for (i,t) in enumerate(data_t)
     integrate(t)
     sim_rv[i] = ps[2].vx
	 chisq += (ps[2].vx-data_rv[i])^2
	 for k in 1:4
	   d_chisq[k] += (ps[2].vx+ps[k*N+2].vx-data_rv[i])^2
	 end
  end
  
  return sim_rv, chisq, d_chisq
end

# Test basic functions
param_init = [4.5,0.6,0.,0.]
setup_rv_simulation(param_init)
sim_rv0 = calc_rvs(param_init)
output = calc_rvs_chisq_grad( param_init )

