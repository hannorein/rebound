/**
 * @file 	problem.c
 * @brief 	Example problem: Restricted three body problem and MPI.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This problem uses MPI to calculate the restricted three
 * body problem. Active particles are copied to all nodes. All other 
 * particles only exist on one node and are not automatically (re-)
 * distributed. There is not domain decomposition used in this example.
 * Run with `mpirun -np 4 nbody`.
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
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "communication_mpi.h"
#include "integrator.h"

int main(int argc, char* argv[]){
	// Setup constants
	integrator		= LEAPFROG;
	boxsize 		= 8; 
	softening		= 1e-6;
	dt 			= 1.0e-2*2.*M_PI;
	N_active 		= 2; 	// Only the star and the planet have non-zero mass
	N_root_x	= 2; N_root_y	= 2; N_root_z	= 1; 
	init_box();
	
	// Initial conditions for star
	struct reb_particle star;
	star.x  = 0; 	star.y  = 0; 	star.z  = 0; 
	star.vx = 0; 	star.vy = 0; 	star.vz = 0;
	star.m  = 1;

	// Initial conditions for planet
	double planet_e = 0.;
	struct reb_particle planet;
	planet.x  = 1.-planet_e; 	planet.y  = 0; 				planet.z  = 0; 
	planet.vx = 0; 			planet.vy = sqrt(2./(1.-planet_e)-1.); 	planet.vz = 0;
	planet.m  = 1e-2;
	
	int _N = 40000;

	// Move to center of mass frame (otherwise planet and star drift out of box)
	double com_x  = (star.x*star.m  + planet.x*planet.m) /(star.m+planet.m);
	double com_y  = (star.y*star.m  + planet.y*planet.m) /(star.m+planet.m);
	double com_z  = (star.z*star.m  + planet.z*planet.m) /(star.m+planet.m);
	double com_vx = (star.vx*star.m + planet.vx*planet.m)/(star.m+planet.m);
	double com_vy = (star.vy*star.m + planet.vy*planet.m)/(star.m+planet.m);
	double com_vz = (star.vz*star.m + planet.vz*planet.m)/(star.m+planet.m);
	planet.x  -= com_x; 	planet.y  -= com_y; 	planet.z  -= com_z;
	planet.vx -= com_vx; 	planet.vy -= com_vy; 	planet.vz -= com_vz;
	star.x    -= com_x; 	star.y    -= com_y; 	star.z    -= com_z;
	star.vx   -= com_vx; 	star.vy   -= com_vy; 	star.vz   -= com_vz;

	// Add active particles on all nodes
	reb_simulation_add(r, star);
	reb_simulation_add(r, planet);
#ifdef MPI
	// Create _N particles in total.
	_N /= mpi_num;
#endif	// MPI
	
	while(N<_N+2){
		double x 	= ((double)rand()/(double)RAND_MAX-0.5)*boxsize*0.9;
		double y 	= ((double)rand()/(double)RAND_MAX-0.5)*boxsize*0.9;
		double a 	= sqrt(x*x+y*y);
		double phi 	= atan2(y,x);
		if (a<.1) continue;
		if (a>boxsize_x/2.*0.9) continue;

		double vkep = sqrt(G*star.m/a);
		struct reb_particle testparticle;
		testparticle.x  = x;
		testparticle.y  = y; 
		testparticle.z  = 1.0e-2*x*((double)rand()/(double)RAND_MAX-0.5);
		testparticle.vx = -vkep*sin(phi);
		testparticle.vy = vkep*cos(phi);
		testparticle.vz = 0;
		testparticle.ax = 0; 
		testparticle.ay = 0; 
		testparticle.az = 0;
		testparticle.m  = 0;
		
		// Add particles locally. This does not distribute particles.
		particles_add_local(testparticle);
	}
}

void heartbeat(struct reb_simulation* r){
	if (reb_simulation_output_check(2.*M_PI)){
		reb_simulation_output_timing();
	}
	if (reb_simulation_output_check(2.*M_PI)){
		reb_simulation_output_ascii("positions.txt");
	}
}

void problem_finish(){
}

