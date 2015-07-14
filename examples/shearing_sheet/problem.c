/**
 * @file 	problem.c
 * @brief 	Example problem: shearing sheet.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example simulates a small patch of Saturn's
 * Rings in shearing sheet coordinates. If you have OpenGL enabled, 
 * you'll see one copy of the computational domain. Press `g` to see
 * the ghost boxes which are used to calculate gravity and collisions.
 * Particle properties resemble those found in Saturn's rings. 
 *
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
#include "rebound.h"
#include "particle.h"
#include "output.h"
#include "tools.h"

double coefficient_of_restitution_bridges(const struct reb_simulation* const r, double v);


void heartbeat(struct reb_simulation* const r);

int main(int argc, char* argv[]) {
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	r->opening_angle2	= .5;					// This determines the precission of the tree code gravity calculation.
	r->integrator			= RB_IT_SEI;
	r->boundary			= RB_BT_SHEAR;
	r->gravity			= RB_GT_TREE;
	r->collision			= RB_CT_TREE;
	double OMEGA 			= 0.00013143527;	// 1/s
	r->ri_sei.OMEGA 		= OMEGA;
	r->G 				= 6.67428e-11;		// N / (1e-5 kg)^2 m^2
	r->softening 			= 0.1;			// m
	r->dt 				= 1e-3*2.*M_PI/OMEGA;	// s
	r->heartbeat			= heartbeat;	// function pointer for heartbeat
	// This example uses two root boxes in the x and y direction. 
	// Although not necessary in this case, it allows for the parallelization using MPI. 
	// See Rein & Liu for a description of what a root box is in this context.
	double surfacedensity 		= 400; 			// kg/m^2
	double particle_density		= 400;			// kg/m^3
	double particle_radius_min 	= 1;			// m
	double particle_radius_max 	= 4;			// m
	double particle_radius_slope 	= -3;	
	double boxsize 			= 100;			// m
	if (argc>1){						// Try to read boxsize from command line
		boxsize = atof(argv[1]);
	}
	reb_configure_box(r, boxsize, 2, 2, 1);
	r->nghostx = 2;
	r->nghosty = 2;
	r->nghostz = 0;
	
	// Initial conditions
	printf("Toomre wavelength: %f\n",4.*M_PI*M_PI*surfacedensity/OMEGA/OMEGA*r->G);
	// Use Bridges et al coefficient of restitution.
	r->coefficient_of_restitution = coefficient_of_restitution_bridges;
	// When two particles collide and the relative velocity is zero, the might sink into each other in the next time step.
	// By adding a small repulsive velocity to each collision, we prevent this from happening.
	r->minimum_collision_velocity = particle_radius_min*OMEGA*0.001;  // small fraction of the shear accross a particle


	// Add all ring paricles
	double total_mass = surfacedensity*r->boxsize.x*r->boxsize.y;
	double mass = 0;
	while(mass<total_mass){
		struct reb_particle pt;
		pt.x 		= reb_tools_uniform(-r->boxsize.x/2.,r->boxsize.x/2.);
		pt.y 		= reb_tools_uniform(-r->boxsize.y/2.,r->boxsize.y/2.);
		pt.z 		= reb_tools_normal(1.);					// m
		pt.vx 		= 0;
		pt.vy 		= -1.5*pt.x*OMEGA;
		pt.vz 		= 0;
		pt.ax 		= 0;
		pt.ay 		= 0;
		pt.az 		= 0;
		double radius 	= reb_tools_powerlaw(particle_radius_min,particle_radius_max,particle_radius_slope);
		pt.r 		= radius;						// m
		double		particle_mass = particle_density*4./3.*M_PI*radius*radius*radius;
		pt.m 		= particle_mass; 	// kg
		reb_add(r, pt);
		mass += particle_mass;
	}
	reb_integrate(r, INFINITY);
}

// This example is using a custom velocity dependend coefficient of restitution
double coefficient_of_restitution_bridges(const struct reb_simulation* const r, double v){
	// assumes v in units of [m/s]
	double eps = 0.32*pow(fabs(v)*100.,-0.234);
	if (eps>1) eps=1;
	if (eps<0) eps=0;
	return eps;
}

void heartbeat(struct reb_simulation* const r){
	if (reb_output_check(r, 1e-3*2.*M_PI/r->ri_sei.OMEGA)){
		reb_output_timing(r, 0);
		//reb_output_append_velocity_dispersion("veldisp.txt");
	}
	if (reb_output_check(r, 2.*M_PI/r->ri_sei.OMEGA)){
		//reb_output_ascii("position.txt");
	}
}

