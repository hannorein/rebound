/**
 * @file 	integrator.c
 * @brief 	Leap-frog integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the leap-frog integration scheme.  
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
#include "particle.h"
#include "main.h"
#include "gravity.h"
#include "boundaries.h"
#include "integrator.h"

// Leapfrog integrator (Drift-Kick-Drift)
// for non-rotating frame.
void integrator_leapfrog_part1(struct Rebound* r){
	const int N = r->N;
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		particles[i].x  += 0.5* r->dt * particles[i].vx;
		particles[i].y  += 0.5* r->dt * particles[i].vy;
		particles[i].z  += 0.5* r->dt * particles[i].vz;
	}
	r->t+=r->dt/2.;
}
void integrator_leapfrog_part2(struct Rebound* r){
	const int N = r->N;
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		particles[i].vx += r->dt * particles[i].ax;
		particles[i].vy += r->dt * particles[i].ay;
		particles[i].vz += r->dt * particles[i].az;
		particles[i].x  += 0.5* r->dt * particles[i].vx;
		particles[i].y  += 0.5* r->dt * particles[i].vy;
		particles[i].z  += 0.5* r->dt * particles[i].vz;
	}
	r->t+=r->dt/2.;
}
	
void integrator_leapfrog_synchronize(struct Rebound* r){
	// Do nothing.
}

void integrator_leapfrog_reset(struct Rebound* r){
	// Do nothing.
}
