/**
 * @brief 	Hybrid symplectic/IAS15 integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Hanno Rein
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "rebound.h"
#include "gravity.h"
#include "integrator.h"
#include "integrator_whfast.h"
#include "integrator_ias15.h"

// Switch to non-symplectic integrator if force_form_star/force_from_other_particle < r->ri_hybrid.switch_ratio.

static double initial_dt = 0;
static unsigned int reb_integrator_hybrid_switch_warning = 0;

static double get_min_ratio(struct reb_simulation* const r){
	const int N = r->N;
	const int N_active = r->N_active;
	const int N_var = r->N_var;
	struct reb_particle* restrict const particles = r->particles;
	struct reb_particle p0 = particles[0];
	const int _N_active = ((N_active==-1)?N:N_active)- N_var;
	const int _N_real   = N - N_var;
	double min_ratio = 1e308;
 

	for (int i=1; i<_N_active; i++){
		struct reb_particle pi = particles[i];
		const double dxi = p0.x - pi.x;
		const double dyi = p0.y - pi.y;
		const double dzi = p0.z - pi.z;
		const double r0i2 = dxi*dxi + dyi*dyi + dzi*dzi;
		const double rhi = r0i2*pow((pi.m/(p0.m*3.)), 2./3.);

		for (int j=1; j<_N_real; j++){
			if (i==j) continue;

			struct reb_particle pj = particles[j];

			const double dx = pi.x - pj.x;
			const double dy = pi.y - pj.y;
			const double dz = pi.z - pj.z;
			const double rij2 = dx*dx + dy*dy + dz*dz;
			const double dxj = p0.x - pj.x;
			const double dyj = p0.y - pj.y;
			const double dzj = p0.z - pj.z;
			const double r0j2 = dxj*dxj + dyj*dyj + dzj*dzj;
			const double rhj = r0j2*pow((pj.m/(p0.m*3.)), 2./3.);

			const double ratio = rij2/(rhi+rhj);

			if (ratio<min_ratio){
				min_ratio = ratio;
			}
		}
	}
	return sqrt(min_ratio);
}


void reb_integrator_hybrid_part1(struct reb_simulation* r){
	const double ratio = get_min_ratio(r);
	if (initial_dt==0.){
		initial_dt = r->dt;
	}
	if (ratio<r->ri_hybrid.switch_ratio){
		if (r->ri_hybrid.mode==SYMPLECTIC){
			reb_integrator_ias15_reset(r); //previous guesses no good anymore
			if (reb_integrator_hybrid_switch_warning==0.){
				reb_integrator_hybrid_switch_warning++;
				fprintf(stderr,"\n\033[1mInfo!\033[0m Switching to HIGHORDER for the first time at t=%.9e.\n",r->t);
			}
			reb_integrator_whfast_synchronize(r);
			r->gravity_ignore_10 = 0;
		}
		r->ri_hybrid.mode = HIGHORDER;
	}else{
		if (r->ri_hybrid.mode==HIGHORDER){
			//reb_integrator_whfast_reset(r); 
			r->ri_whfast.recalculate_jacobi_this_timestep = 1;
			r->dt = initial_dt;
		}
		r->ri_hybrid.mode = SYMPLECTIC;
	}
	switch(r->ri_hybrid.mode){
		case SYMPLECTIC:
			reb_integrator_whfast_part1(r);
			break;
		case HIGHORDER:
			reb_integrator_ias15_part1(r);
			break;
		default:
			break;
	}
}
void reb_integrator_hybrid_part2(struct reb_simulation* r){
	switch(r->ri_hybrid.mode){
		case SYMPLECTIC:
			reb_integrator_whfast_part2(r);
			break;
		case HIGHORDER:
			reb_integrator_ias15_part2(r);
			break;
		default:
			break;
	}
}
	
void reb_integrator_hybrid_synchronize(struct reb_simulation* r){
	switch(r->ri_hybrid.mode){
		case SYMPLECTIC:
			reb_integrator_whfast_synchronize(r);
			break;
		default:
			break;
	}
}

void reb_integrator_hybrid_reset(struct reb_simulation* r){
	r->ri_hybrid.mode = SYMPLECTIC;
	reb_integrator_hybrid_switch_warning = 0;
	reb_integrator_whfast_reset(r);
	reb_integrator_ias15_reset(r);
	r->ri_hybrid.switch_ratio = 8.;
	initial_dt = 0.;
}
