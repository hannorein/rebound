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

extern double OMEGA;

#ifndef COLLISIONS_NONE
extern double coefficient_of_restitution;
extern double minimum_collision_velocity;
extern double (*coefficient_of_restitution_for_velocity)(double); 
double coefficient_of_restitution_bridges(double v){
	// assumes v in units of [m/s]
	double eps = 0.32*pow(fabs(v)*100.,-0.234);
	if (eps>1) eps=1;
	if (eps<0) eps=0;
	return eps;
}
#endif // COLLISIONS_NONE
#ifdef GRAVITY_GRAPE
extern double gravity_range;
#endif // GRAVITY_GRAPE


void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA 				= 0.00013143527;	// 1/s
	G 				= 6.67428e-11;		// N / (1e-5 kg)^2 m^2
	dt 				= 1e-3*2.*M_PI/OMEGA;	// s
	root_nx = 10; root_ny = 1; root_nz = 1;
	nghostx = 1; nghosty = 1; nghostz = 0; 			// Use two one ring (+cutoff, see below)
	double surfacedensity 		= 400; 			// kg/m^2
	double particle_density		= 400;			// kg/m^3
	double particle_radius_min 	= 1;			// m
	double particle_radius_max 	= 4;			// m
	double particle_radius_slope 	= -3;	
	boxsize 			= 100;
	if (argc>1){						// Try to read boxsize from command line
		boxsize = atof(argv[1]);
	}
	init_box();
#ifdef GRAVITY_GRAPE
	gravity_range = boxsize/2.;
#endif // GRAVITY_GRAPE
	
	// Initial conditions
	printf("Toomre wavelength: %f\n",2.*M_PI*M_PI*surfacedensity/OMEGA/OMEGA*G);
#ifndef COLLISIONS_NONE
	// Use Bridges et al coefficient of restitution.
	coefficient_of_restitution_for_velocity	= coefficient_of_restitution_bridges;
	minimum_collision_velocity		= particle_radius_min*OMEGA*0.001;  // small fraction of the shear
	softening 				= 0.1;			// m
#else  // COLLISIONS_NONE
	softening				= particle_radius_max;
#endif // COLLISIONS_NONE
	double total_mass = surfacedensity*boxsize_x*boxsize_y;
	double mass = 0;
	while(mass<total_mass){
		struct particle pt;
		pt.x 		= tools_uniform(-boxsize_x/2.,boxsize_x/2.);
		pt.y 		= tools_uniform(-boxsize_y/2.,boxsize_y/2.);
		pt.z 		= tools_normal(1.);					// m
		pt.vx 		= 0;
		pt.vy 		= -1.5*pt.x*OMEGA;
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
}

void problem_inloop(){
}

void output_ascii_mod(char* filename){
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,mpi_id);
	FILE* of = fopen(filename_mpi,"w"); 
#else // MPI
	FILE* of = fopen(filename,"w"); 
#endif // MPI
	if (of==NULL){
		printf("\n\nError while opening file '%s'.\n",filename);
		return;
	}
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		fprintf(of,"%e\t%e\t%e\t%e\n",p.x,p.y,p.z,p.r);
	}
	fclose(of);
}

void problem_output(){
	if (output_check(1e-1*2.*M_PI/OMEGA)){
		output_timing();
	}
	if (output_check(.2*M_PI/OMEGA)){
		output_ascii_mod("ascii.txt");
	}
}

void problem_finish(){
}
