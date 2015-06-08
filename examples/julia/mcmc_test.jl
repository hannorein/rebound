if !isdefined(:LIB_REBOUND_INITIALIZED)
  include("librebound.jl")
  using Rebound
  LIB_REBOUND_INITIALIZED = true
end

using Distributions
using Lora
# using PyPlot

# Setting up simulations
function setup_rv_simulation( param::Array{Float64,1} )
  @assert(length(param)==4)
  star_mass = 1.0
  a=param[1]
  ecc=param[2]
  omega=param[3]
  anom=param[4]
  m_pl = 1e-3
  Rebound.reset()
  #set_integrator("whfast")
  set_dt(0.1)        # in year/(2*pi)

  Rebound.add( rebound_particle_basic(1.0) )
  Rebound.add( Rebound.tools_init_orbit2d(star_mass, m_pl, a, ecc, omega, anom) )
end

# Calc RVs
function calc_rvs( data_t::Array{Float64,1} )
  @assert(length(data_t)>=1)
  ps = Rebound.get_particles()
  sim_rv = zeros(length(data_t))
  for (i,t) in enumerate(data_t)
     integrate(t)
     sim_rv[i] = ps[2].vx
  end
  return sim_rv
end

# Calc RVs, chisq and its gradient
function calc_rvs_chisq_grad( param::Array{Float64,1}, data_t::Vector{Float64}, data_rv::Vector{Float64}  )
  @assert(length(param)==4)
  @assert(length(data_t)==length(data_rv)>=3)

  star_mass = 1.0
  # a=param[1]; ecc=param[2]; omega=param[3]; anom=param[4]
  m_pl = 1e-3
  Rebound.reset()
  Rebound.set_N_megnopp(4)
  #set_integrator("whfast")
  set_dt(0.1)        # in year/(2*pi)

  Rebound.add( rebound_particle_basic(star_mass) )
  Rebound.add( Rebound.tools_init_orbit2d(star_mass, m_pl, param[1], param[2], param[3], param[4]) )
  Rebound.tools_megno_init(0.)
  ps = Rebound.get_particles()
  num_real_part = 2
  for k in [1:4]
    delta = 1e-6
	param2 = param
	param2[k] += delta

	mp = Rebound.tools_init_orbit2d(star_mass, m_pl, param2[1], param2[2], param2[3], param2[4]) 
	vari = (num_real_part*k)+2
	vp = rebound_particle_basic(ps[2].x - mp.x, ps[2].y - mp.y, ps[2].z - mp.z,  
			       ps[2].vx - mp.vx, ps[2].vy - mp.vy, ps[2].vz - mp.vz, 
			       0.0, 0.0, 0.0, m_pl)
	# println("# N= ",get_N(), " vari= ",vari)
        ps[vari] = vp
  end

  sim_rv = zeros(length(data_t))
  chisq = 0.0
  d_chisq = zeros(4) # length(param))
  for (i,t) in enumerate(data_t)
     integrate(t)
     sim_rv[i] = ps[2].vx 
     chisq += (ps[2].vx-data_rv[i])^2
     for k in 1:4
       d_chisq[k] += (ps[2].vx+ps[k*num_real_part+2].vx-data_rv[i])^2
     end
  end
 
  return sim_rv, chisq, d_chisq
end

function logl( param::Vector{Float64} )
  -0.5*calc_rvs_chisq_grad( param, data_t, data_rv)[2]
end

function grad_logl( param::Vector{Float64} )
  -0.5*calc_rvs_chisq_grad( param, data_t, data_rv)[3]
end

function logl_and_grad_logl( param::Vector{Float64} )
  rvs, chisq, d_chisq = calc_rvs_chisq_grad( param, data_t, data_rv)
  loglikelihood = -0.5*chisq
  grad_ll = -0.5.*d_chisq
  return loglikelihood, grad_ll
end


# Simulated Data:
data_N = 30
data_t = linspace(0.1,100,data_N).+0.5.*rand(data_N)
# data_rv = 0.500883379885605.*(sin(2pi.*data_t/50.+0.2*pi).+0.6.*(rand(data_N).-0.5))
# plot(data_t,data_rv)

param_true = [4.0,0.6,0.5*pi,pi]
setup_rv_simulation(param_true)
println("# Calculating rvs")
tic()
data_true = calc_rvs(data_t)
data_rv = data_true .+ 0.001*randn(length(data_true))
toc()
println("# Calculating: output = calc_rvs_chisq_grad( param_init ) ")
tic()
output_at_param_true = calc_rvs_chisq_grad( param_true, data_t, data_rv )
toc()

param_init = copy(param_true)
param_init = [4.1,0.5,0.4*pi,0.9*pi]
println("# Constructing statistical model")
tic()
mcmodel = model(logl, grad=grad_logl, init=param_init )
toc()
println("# Running HMC")
tic()
mcchain = run(mcmodel, HMC(10,0.1), SerialMC(1:1000))
toc()
println("# Acceptance rate = ", acceptance(mcchain))
describe(mcchain)


