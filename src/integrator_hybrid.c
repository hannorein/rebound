/**
 * @file 	integrator.c
 * @brief 	Hybrid symplectic/IAS15 integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the leap-frog integration scheme.  
 * This scheme is second order accurate, symplectic and well suited for 
 * non-rotating coordinate systems. Note that the scheme is formally only
 * first order accurate when velocity dependent forces are present.
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

// Switch to non-symplectic integrator if force_form_star/force_from_other_particle < reb_integrator_hybrid_switch_ratio.
// Default corresponds to about 10 Hill Radii 
double reb_integrator_hybrid_switch_ratio = 100; 

static double initial_dt = 0;
static unsigned int reb_integrator_hybrid_switch_warning = 0;

static double get_min_ratio(struct reb_context* const r){
	const int N = r->N;
	const int N_active = r->N_active;
	const int N_megno = r->N_megno;
	struct reb_particle* restrict const particles = r->particles;
	struct reb_particle p0 = particles[0];
	const int _N_active = ((N_active==-1)?N:N_active)- N_megno;
	const int _N_real   = N - N_megno;
	double min_ratio = 1e308;
	for (int i=1; i<_N_active; i++){
		struct reb_particle pi = particles[i];
	for (int j=1; j<_N_real; j++){
		if (i==j) continue;
		const double dxj = p0.x - particles[j].x;
		const double dyj = p0.y - particles[j].y;
		const double dzj = p0.z - particles[j].z;
		const double r0j2 = dxj*dxj + dyj*dyj + dzj*dzj;

		const double dx = pi.x - particles[j].x;
		const double dy = pi.y - particles[j].y;
		const double dz = pi.z - particles[j].z;
		const double rij2 = dx*dx + dy*dy + dz*dz;
		
		const double F0j = p0.m/r0j2;
		const double Fij = pi.m/rij2;

		const double ratio = F0j/Fij;
			
		if (ratio<min_ratio){
			min_ratio = ratio;
		}
	}
	}
	return min_ratio;
}


static double ratio;
enum {SYMPLECTIC, HIGHORDER} reb_integrator_hybrid_mode = SYMPLECTIC;
void reb_integrator_hybrid_part1(struct reb_context* r){
	ratio = get_min_ratio(r);
	if (initial_dt==0.){
		initial_dt = r->dt;
	}
	if (ratio<reb_integrator_hybrid_switch_ratio){
		if (reb_integrator_hybrid_mode==SYMPLECTIC){
			reb_integrator_ias15_reset(r); //previous guesses no good anymore
			if (reb_integrator_hybrid_switch_warning==0.){
				reb_integrator_hybrid_switch_warning++;
				fprintf(stderr,"\n\033[1mInfo!\033[0m Switching to HIGHORDER for the first time at t=%.9e.\n",r->t);
			}
			reb_integrator_whfast_synchronize(r);
			r->gravity_ignore_10 = 0;
		}
		reb_integrator_hybrid_mode = HIGHORDER;
	}else{
		if (reb_integrator_hybrid_mode==HIGHORDER){
			//reb_integrator_whfast_reset(r); 
			r->ri_whfast.recalculate_jacobi_this_timestep = 1;
			r->dt = initial_dt;
		}
		reb_integrator_hybrid_mode = SYMPLECTIC;
	}
	switch(reb_integrator_hybrid_mode){
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
void reb_integrator_hybrid_part2(struct reb_context* r){
	switch(reb_integrator_hybrid_mode){
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
	
void reb_integrator_hybrid_synchronize(struct reb_context* r){
	switch(reb_integrator_hybrid_mode){
		case SYMPLECTIC:
			reb_integrator_whfast_synchronize(r);
			break;
		default:
			break;
	}
}

void reb_integrator_hybrid_reset(struct reb_context* r){
	reb_integrator_hybrid_mode = SYMPLECTIC;
	reb_integrator_hybrid_switch_warning = 0;
	reb_integrator_whfast_reset(r);
	reb_integrator_ias15_reset(r);
	reb_integrator_hybrid_switch_ratio = 400.;
	initial_dt = 0.;
}
