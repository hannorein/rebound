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

// Slightly dirty trick to rename function for librebound use
#ifdef LIBREBOUND
#define integrator_part1                      integrator_leapfrog_part1
#define integrator_part2                      integrator_leapfrog_part2
#define integrator_force_is_velocitydependent integrator_leapfrog_force_is_velocitydependent
#define integrator_epsilon                    integrator_leapfrog_epsilon
#define integrator_min_dt                     integrator_leapfrog_min_dt
#define integrator_reset                      integrator_leapfrog_reset
#endif // LIBREBOUND

// These variables have no effect for leapfrog.
int integrator_force_is_velocitydependent 	= 1;
double integrator_epsilon 			= 0;
double integrator_min_dt 			= 0;


// Leapfrog integrator (Drift-Kick-Drift)
// for non-rotating frame.
void integrator_part1(){
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		particles[i].x  += 0.5* dt * particles[i].vx;
		particles[i].y  += 0.5* dt * particles[i].vy;
		particles[i].z  += 0.5* dt * particles[i].vz;
	}
	t+=dt/2.;
}
void integrator_part2(){
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		particles[i].vx += dt * particles[i].ax;
		particles[i].vy += dt * particles[i].ay;
		particles[i].vz += dt * particles[i].az;
		particles[i].x  += 0.5* dt * particles[i].vx;
		particles[i].y  += 0.5* dt * particles[i].vy;
		particles[i].z  += 0.5* dt * particles[i].vz;
	}
	t+=dt/2.;
}
	

void integrator_reset(){
	// Do nothing.
}
