/**
 * @file 	problem.c
 * @brief 	Example problem: shearing sheet.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This problem is identical to the other shearing
 * sheet examples but uses an FFT based gravity solver. 
 * To run this example, you need to install the FFTW library. 
 * Collisions are detected using a plane sweep algorithm. 
 * There is no tree present in this simulation.
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
#include <sys/time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "communication_mpi.h"
#include "tree.h"
#include "tools.h"
#include "integrator.h"

extern double OMEGA;
extern double OMEGAZ;
extern double coefficient_of_restitution;
extern double minimum_collision_velocity;

extern double (*coefficient_of_restitution_for_velocity)(double); 
double coefficient_of_restitution_bridges(double v); 

int main(int argc, char* argv[]){
	// Setup constants
	integrator			= SEI;
	OMEGA 				= 0.00013143527;	// 1/s
	OMEGAZ 				= 3.6*0.00013143527;	// 1/s
	G 				= 6.67428e-11;		// N / (1e-5 kg)^2 m^2
	softening 			= 0.1;			// m
	dt 				= 1e-3*2.*M_PI/OMEGA;	// s
	int ngrid 			= 64;
	N_root_x = ngrid; N_root_y = ngrid; N_root_z = ngrid/2;
	double surfacedensity 		= 400; 			// kg/m^2
	double particle_density		= 400;			// kg/m^3
	double particle_radius_min 	= 1;			// m
	double particle_radius_max 	= 4;			// m
	double particle_radius_slope 	= -3;	
	boxsize 			= 200/(double)ngrid;
	if (argc>1){						// Try to read boxsize from command line
		boxsize = atof(argv[1]);
	}
	init_box();
	
	// Use Bridges et al coefficient of restitution.
	coefficient_of_restitution_for_velocity = coefficient_of_restitution_bridges;
	minimum_collision_velocity = particle_radius_min*OMEGA*0.001;  // small fraction of the shear
	double total_mass = surfacedensity*boxsize_x*boxsize_y;
#ifdef MPI
	// Only initialise particles on master. This should also be parallelied but the details depend on the individual problem.
	if (mpi_id==0){
#endif
		double mass = 0;
		while(mass<total_mass){
			struct reb_particle pt;
			pt.x 		= reb_random_uniform(-boxsize_x/2.,boxsize_x/2.);
			pt.y 		= reb_random_uniform(-boxsize_y/2.,boxsize_y/2.);
			pt.z 		= reb_random_normal(1.);					// m
			pt.vx 		= 0;
			pt.vy 		= -1.5*pt.x*OMEGA;
			pt.vz 		= 0;
			pt.ax 		= 0;
			pt.ay 		= 0;
			pt.az 		= 0;
			double radius 	= reb_random_powerlaw(particle_radius_min,particle_radius_max,particle_radius_slope);
#ifndef COLLISIONS_NONE
			pt.r 		= radius;						// m
#endif
			double		particle_mass = particle_density*4./3.*M_PI*radius*radius*radius;
			pt.m 		= particle_mass; 	// kg
			reb_simulation_add(r, pt);
			mass += particle_mass;
		}
#ifdef MPI
	}
#endif
}

double coefficient_of_restitution_bridges(double v){
	// assumes v in units of [m/s]
	double eps = 0.32*pow(fabs(v)*100.,-0.234);
	if (eps>1) eps=1;
	if (eps<0) eps=0;
	return eps;
}

void heartbeat(struct reb_simulation* r){
	if (reb_simulation_output_check(10.0*dt)){
		reb_simulation_output_timing();
	}
}

void problem_finish(){
}
