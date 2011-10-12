/**
 * @file 	problem.c
 * @brief 	Example problem: shearing sheet.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This problem uses shearing sheet boundary
 * conditions. Particle properties resemble those found in 
 * Saturn's rings. 
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

#ifdef INTEGRATOR_SEI
extern double OMEGA;
extern double OMEGAZ;
#endif 	// INTEGRATOR_SEI
extern double coefficient_of_restitution;
extern double minimum_collision_velocity;

extern double (*coefficient_of_restitution_for_velocity)(double); 
double coefficient_of_restitution_bridges(double v); 

void problem_init(int argc, char* argv[]){
	// Setup constants
#ifdef INTEGRATOR_SEI
	OMEGA 				= 0.00013143527;	// 1/s
	OMEGAZ 				= 3.6*0.00013143527;	// 1/s
#endif 	// INTEGRATOR_SEI
	G 				= 6.67428e-11;		// N / (1e-5 kg)^2 m^2
	softening 			= 0.1;			// m
#ifdef INTEGRATOR_SEI
	dt 				= 1e-3*2.*M_PI/OMEGA;	// s
#else 	// INTEGRATOR_SEI
	dt 				= 1e2;			// s
#endif 	// INTEGRATOR_SEI
	int ngrid 			= 64;
	root_nx = ngrid; root_ny = ngrid; root_nz = ngrid/2;
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
#ifdef INTEGRATOR_SEI
	minimum_collision_velocity = particle_radius_min*OMEGA*0.001;  // small fraction of the shear
#else	// INTEGRATOR_SEI
	minimum_collision_velocity = particle_radius_min/dt*0.001;  // small fraction of the shear
#endif	// INTEGRATOR_SEI
	double total_mass = surfacedensity*boxsize_x*boxsize_y;
#ifdef MPI
	// Only initialise particles on master. This should also be parallelied but the details depend on the individual problem.
	if (mpi_id==0){
#endif
		double mass = 0;
		while(mass<total_mass){
			struct particle pt;
			pt.x 		= tools_uniform(-boxsize_x/2.,boxsize_x/2.);
			pt.y 		= tools_uniform(-boxsize_y/2.,boxsize_y/2.);
			pt.z 		= tools_normal(1.);					// m
			pt.vx 		= 0;
#ifdef INTEGRATOR_SEI
			pt.vy 		= -1.5*pt.x*OMEGA;
#else	// INTEGRATOR_SEI
			pt.vy 		= 0;
#endif 	// INTEGRATOR_SEI
			pt.vz 		= 0;
			pt.ax 		= 0;
			pt.ay 		= 0;
			pt.az 		= 0;
			double radius 	= tools_powerlaw(particle_radius_min,particle_radius_max,particle_radius_slope);
#ifndef COLLISIONS_NONE
			pt.r 		= radius;						// m
#endif
			double		particle_mass = particle_density*4./3.*M_PI*radius*radius*radius;
			pt.m 		= particle_mass; 	// kg
			particles_add(pt);
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

void problem_inloop(){
}

void problem_output(){
	if (output_check(10.0*dt)){
		output_timing();
	}
}

void problem_finish(){
}
