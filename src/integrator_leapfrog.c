/**
 * @file 	integrator.c
 * @brief 	Leap-frog integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	This file implements the leap-frog integration scheme.  
 * This scheme is second order accurate, symplectic and well suited for 
 * non-rotating coordinate systems. Note that the scheme is formally only
 * first order accurate when velocity dependent forces are present.
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
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
#include "rebound.h"

// Leapfrog integrator (Drift-Kick-Drift)
// for non-rotating frame.
void reb_integrator_leapfrog_part1(struct reb_simulation* r){
    r->gravity_ignore_terms = 0;
	const int N = r->N;
	struct reb_particle* restrict const particles = r->particles;
	const double dt = r->dt;
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		particles[i].x  += 0.5* dt * particles[i].vx;
		particles[i].y  += 0.5* dt * particles[i].vy;
		particles[i].z  += 0.5* dt * particles[i].vz;
	}
	r->t+=dt/2.;
}
void reb_integrator_leapfrog_part2(struct reb_simulation* r){
	const int N = r->N;
	struct reb_particle* restrict const particles = r->particles;
	const double dt = r->dt;
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		particles[i].vx += dt * particles[i].ax;
		particles[i].vy += dt * particles[i].ay;
		particles[i].vz += dt * particles[i].az;
		particles[i].x  += 0.5* dt * particles[i].vx;
		particles[i].y  += 0.5* dt * particles[i].vy;
		particles[i].z  += 0.5* dt * particles[i].vz;
	}
	r->t+=dt/2.;
	r->dt_last_done = r->dt;
}
	
void reb_integrator_leapfrog_synchronize(struct reb_simulation* r){
	// Do nothing.
}

void reb_integrator_leapfrog_reset(struct reb_simulation* r){
	// Do nothing.
}
