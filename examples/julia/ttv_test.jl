if !isdefined(:LIB_REBOUND_INITIALIZED)
  include("librebound.jl")
  using Rebound
  LIB_REBOUND_INITIALIZED = true
end

function est_ttv_newton()
  local part = Rebound.get_particles()
  local dx = [ part[2].x-part[1].x,   part[2].y-part[1].y,   part[2].z-part[1].z  ]
  local dv = [ part[2].vx-part[1].vx, part[2].vy-part[1].vy, part[2].vz-part[1].vz ]
  local da = [ part[2].ax-part[1].ax, part[2].ay-part[1].ay, part[2].az-part[1].az ]
  local dt
  if false # general
   const los = [ 0.0, 0.0, 1.0 ]
   local dxcl = cross(dx,los)
   local dvcl = cross(dv,los)
   local dacl = cross(da,los)
   #b2 = dot(dxcl,dxcl)
   local db2dt = 2*dot(dxcl,dvcl)
   local d2b2dt = 2*(dot(dxcl,dacl)+dot(dvcl,dvcl))
   dt = -db2dt/d2b2dt
  else # los on z-axis
   local D = dx[1]*dv[1]+dx[2]*dv[2]
   local Ddot = dv[1]^2+dv[2]^2+dx[1]*da[1]+dx[2]*da[2]
   dt = -D/Ddot
   #b2 = (dx[1]+dt*(dv[1]+0.5*dt*da[1]))^2 + (dx[2]+dt*(dv[2]+0.5*dt*da[2]))^2
  end
  return dt
end

function calc_ttvs_rebound(log_times::Vector{Float64} )
	@assert(length(log_times)>=3)
	local ttvs = zeros(length(log_times))
	local tts = zeros(length(log_times))

    	# local e_init = Rebound.tools_energy()
        local time_step = copysign(get_dt(),log_times[1]-get_t())
 	tic()
	for (i,t) in enumerate(log_times)
	   # println("# ",i," t = ",t)
	   local tguess = t
	   local  flag = integrate(tguess,1) # exactFinishTime=1 for exact end times
	     @assert(flag==0)
 	     local dt1 = est_ttv_newton()
	     local tt1 = get_t() + dt1
	     # println("# i=$i dt=", dt1, "  tt=",tt1)
             time_step = copysign(time_step,dt1)
	     set_dt(time_step)
	     tguess = tt1

	     flag = integrate(tguess,1) # exactFinishTime=1 for exact end times
	     @assert(flag==0)
 	     local dt2 = est_ttv_newton()
             local tt2 = get_t() + dt2
	     # println("# i=$i dt=", dt2, "  tt=",tt2)

	    tts[i] = (dt1*tt2-dt2*tt1)/(dt1-dt2)
	    ttvs[i] = tts[i]-t
 	   time_step = abs(time_step)
	   set_dt(time_step)
	end
	toc()
    	# local e_final = Rebound.tools_energy()
    	# local delta_e = e_final - e_init
    	# println("# t= ", get_t(), " Integrated for ", last(tts)-first(tts), " dE/E = ",delta_e/e_init)
	return tts, ttvs
end

# For system based on KOI-142, but with m3 divided by 10
function setup_test_calc_ttvs_rebound()
        Rebound.reset()
	set_integrator("whfast")
        # Initial conditions from TTVFast's example
        G = 0.000295994511  # to match TTVfast units of days, AU, Msol
	set_G(G)
	star_mass = 0.95573417954 
        pstar = rebound_particle_basic(star_mass) #  Star
        Rebound.add(pstar);  
        Rebound.add(rebound_particle_basic(4.2751105789149389e-03,-1.5242519870492784e-03,9.4799180429814917e-02,-5.4584946596113952e-02,9.7656270156749417e-06,-6.0736246062660626e-04,m=0.00002878248) )
        #Rebound.add(rebound_particle_basic(1.3559014131822117e-01,-1.0049182503274595e-03,-5.0033242877868256e-02,1.4859847408958543e-02,1.9180380744132197e-03,4.2870425084428648e-02,m=0.00061895914) ) # with actual masses, ttv's are too big for this algorithm
        Rebound.add(rebound_particle_basic(1.3559014131822117e-01,-1.0049182503274595e-03,-5.0033242877868256e-02,1.4859847408958543e-02,1.9180380744132197e-03,4.2870425084428648e-02,m=0.000061895914) )

	t_init = -1045.0
	set_t(t_init)
	dt = 0.054 # in days, not years or year/2pi
	dt = 0.1 # in days, not years or year/2pi
	set_dt(dt)
#        tr1_first = -1.044921707263204e+03  # for real masses
# 	tr1_last  = 1.692353532433625e+03
#        num_tr1 = 251
        tr1_first = -1.044921707507387e+03  # for m_3 -> m_3/10
        tr1_last = 1.696433128728580e+03
        num_tr1 = 252
        
	Period = (tr1_last-tr1_last)/(num_tr1-1)
	#println("# t_init=$t_init P=$Period t_first=$tr1_first")
	log_times = linspace(tr1_first,tr1_last,num_tr1)
end


tt_guesses = setup_test_calc_ttvs_rebound()
ttv_output = calc_ttvs_rebound(tt_guesses)
writedlm("koi142_rebound.out", ttv_output);

