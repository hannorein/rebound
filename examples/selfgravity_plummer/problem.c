/**
 * @file 	problem.c
 * @brief 	Example problem: self-gravity disc.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	A self-gravitating plummer sphere is integrated using
 * the leap frog integrator. Collisions are not resolved.
 * 
 * @section 	LICENSE
 * Copyright (c) 2013 Hanno Rei
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
#include <sys/time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "communication_mpi.h"
#include "tree.h"
#include "tools.h"

extern double opening_angle2;
extern int Nmax;

void problem_init(int argc, char* argv[]){
	// Setup constants
	G 		= 1;		

	// Setup particles
	int _N = 100; 			// Number of particles
	double M = 1;			// Total mass of the cluster
	double R = 1;			// Radius of the cluster
	double E = 3./64.*M_PI*M*M/R;	// Energy of the cluster
	double r0 = 16./(3.*M_PI)*R;	// Chacateristic length scale
	double t0 = G*pow(M,5./2.)*pow(4.*E,-3./2.)*(double)_N/log(0.4*(double)_N); // Rellaxation time
	printf("Characteristic size:              %f\n", r0);
	printf("Characteristic time (relaxation): %f\n", t0);

	dt 		= 1e-4*t0; 	// timestep
	tmax		= 20.*t0;	// Integration time
	boxsize 	= 20.*r0;	// Size of box (open boundaries)
	softening 	= 0.01*r0;	// Softening parameter

	init_box();			// Set up box (uses variable 'boxsize')
	
	tools_init_plummer(_N, M, R);	// Add particles
	
	tools_move_to_center_of_momentum(); // Move to rest frame 
}

void problem_inloop(){
}

void problem_output(){
	if (output_check(10.0*dt)) output_timing();
}

void problem_finish(){
}
