/**
 * @file 	integrator.c
 * @brief 	RADAU15 integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the radau15 integration scheme.  
 * This scheme is a fifteenth order integrator well suited for 
 * high accuracy orbit integration with non-conservative forces.
 * See Everhart, 1985, ASSL Vol. 115, IAU Colloq. 83, Dynamics of 
 * Comets, Their Origin and Evolution, 185.
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Dave Spiegel.
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
#include "problem.h"
#ifdef TREE
#error RADAU15 integrator not working with TREE module.
#endif
#ifdef MPI
#error RADAU15 integrator not working with MPI.
#endif


void integrator_part1(){
	// Do nothing here. This is for the first drift part in a leapfrog-like DKD integrator.
}

// This function updates the acceleration on all particles. 
// It uses the current position and velocity data in the (struct particle*) particles structure.
// Note: this does currently not work with MPI or any TREE module.
void integrator_update_acceleration(){
	gravity_calculate_acceleration();
	if (problem_additional_forces) problem_additional_forces();
}

void integrator_part2(){
	/////////////////////////////////////////////////////////////////////////////////////////////////
	// Available variables at the beginning of this function:
	// t: 			current time
	// dt: 			suggested times step (doesn't have to be fixed)
	// particles:		x,y,z,vx,vy,vz,ax,ay,az of all particles at the beginning of the timestep.
	/////////////////////////////////////////////////////////////////////////////////////////////////




	// Do the real work here.
	
	// Note that when calling integrator_update_acceleration() within this function, the 
	// correct position and velocities must be stored in the particles[] array. 
	// The updated accelerations will then also be stored in the particles[] array.
	// It's like if the function F(y,y',t) can only take particles[] as an argument, nothing else.
	// You may, however, change the pointer to the particles array temporarily to help you with that. 
	// For example do something like this:
	struct particle* old_particles = particles;		// Save old particle array pointer.
	particles = malloc(sizeof(struct particle)*N);		// Create space for temporary array.
	for (int i=0;i<N;i++){
		particles[i] = old_particles[i];		// Copy old particle data.

		// Let's do an implicit backward Euler method as an example.
		// Note: this is a bad integrator which need a really small timestep.
		// It's not time-reversible and particles fall into the star quickly.
		particles[i].x  += dt*particles[i].vx;		// First guess.	
		particles[i].y  += dt*particles[i].vy;		
		particles[i].z  += dt*particles[i].vz;		
		particles[i].vx += dt*particles[i].ax;		
		particles[i].vy += dt*particles[i].ay;		
		particles[i].vz += dt*particles[i].az;		
	}
	for (int iterations=0;iterations<1000;iterations++){
		integrator_update_acceleration(); // This updates ax,ay,az in the particles[] array.
		for (int i=0;i<N;i++){
			particles[i].x  = old_particles[i].x  + dt*particles[i].vx;	// Each iter	
			particles[i].y  = old_particles[i].y  + dt*particles[i].vy;		
			particles[i].z  = old_particles[i].z  + dt*particles[i].vz;		
			particles[i].vx = old_particles[i].vx + dt*particles[i].ax;		
			particles[i].vy = old_particles[i].vy + dt*particles[i].ay;		
			particles[i].vz = old_particles[i].vz + dt*particles[i].az;		
		}
	}
	for (int i=0;i<N;i++){
		old_particles[i] = particles[i];		// Copy particles back to the old particle data structure.
	}
	free(particles);					// Release memory of temporary array.
	particles = old_particles;				// Change pointer back to old particle array pointer.


	// Advance the time. Because the function is using an adaptive timestep, use the dt that has been computed self-consistely.
	t+=dt; 	



	/////////////////////////////////////////////////////////////////////////////////////////////////
	// These vailables should be set at the end of the timestep:
	// dt: 			actual timestep done
	// particles:		x,y,z,vx,vy,vz of all particles at the end of the timestep.
	/////////////////////////////////////////////////////////////////////////////////////////////////
}
	

