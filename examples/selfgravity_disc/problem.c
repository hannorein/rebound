/**
 * @file 	problem.c
 * @brief 	Example problem: self-gravity disc.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	A self-gravitating disc is integrated using
 * the leap frog integrator. Collisions are not resolved.
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
#include "tools.h"
#include "output.h"


void heartbeat(struct reb_simulation* const r);

int main(int argc, char* argv[]){
	struct reb_simulation* const r = reb_create_simulation();
	// Setup constants
	r->integrator	= RB_IT_LEAPFROG;
	r->gravity	= RB_GT_TREE;
	r->boundary	= RB_BT_OPEN;
	r->opening_angle2	= 1.5;		// This constant determines the accuracy of the tree code gravity estimate.
	r->G 		= 1;		
	r->softening 	= 0.02;		// Gravitational softening length
	r->dt 		= 3e-2;		// Timestep
	const double boxsize = 10.2;
	reb_configure_box(r,boxsize,1,1,1);

	// Setup particles
	double disc_mass = 2e-1;	// Total disc mass
	int N = 10000;			// Number of particles
	// Initial conditions
	struct reb_particle star;
	star.x 		= 0; star.y 	= 0; star.z	= 0;
	star.vx 	= 0; star.vy 	= 0; star.vz 	= 0;
	star.m 		= 1;
	reb_add(r, star);
	for (int i=0;i<N;i++){
		struct reb_particle pt;
		double a	= reb_tools_powerlaw(boxsize/10.,boxsize/2./1.2,-1.5);
		double phi 	= reb_tools_uniform(0,2.*M_PI);
		pt.x 		= a*cos(phi);
		pt.y 		= a*sin(phi);
		pt.z 		= a*reb_tools_normal(0.001);
		double mu 	= star.m + disc_mass * (pow(a,-3./2.)-pow(boxsize/10.,-3./2.))/(pow(boxsize/2./1.2,-3./2.)-pow(boxsize/10.,-3./2.));
		double vkep 	= sqrt(r->G*mu/a);
		pt.vx 		=  vkep * sin(phi);
		pt.vy 		= -vkep * cos(phi);
		pt.vz 		= 0;
		pt.m 		= disc_mass/(double)N;
		reb_add(r, pt);
	}

	r->heartbeat = heartbeat;
	reb_integrate(r,0);
}

void heartbeat(struct reb_simulation* const r){
	if (reb_output_check(r,10.0*r->dt)) reb_output_timing(r,0);
}
