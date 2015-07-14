/**
 * @file 	problem.c
 * @brief 	Example problem: restarting simulation.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example demonstrates how to restart a simulation
 * using a binary file. A shearing sheet ring simulation is used, but
 * the same method can be applied to any other type of simulation.
 * 
 * First, run the program with `./rebound`.
 * Random initial conditions are created and
 * a restart file is written once per orbit.
 * Then, to restart the simulation, run the 
 * program with `./rebound --restart restart.bin`.
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
#include "rebound.h"
#include "output.h"
#include "input.h"

void heartbeat(struct Rebound* const r);

int main(int argc, char* argv[]){
	// Setup constants
	struct Rebound* r;
	
	const int restart = 1;

	if (restart){
		r = rebound_init_from_binary("restart.bin");
	}else{
		r = rebound_init();

		r->integrator	= SEI;
		r->collision	= RB_CT_DIRECT;
		r->ri_sei.OMEGA	= 1.;	
		r->dt 		= 1e-4*2.*M_PI; 
		rebound_configure_box(r,1.,1,1,1);
		r->coefficient_of_restitution = 0.5;
		r->minimum_collision_velocity = 0.01;
		r->nghostx = 1; 
		r->nghosty = 1; 
		r->nghostz = 0;

		// Check if simulation is restarted
		// Initial conditions
		while (r->N<50){
			struct Particle p;
			p.x  = ((double)rand()/(double)RAND_MAX-0.5)*r->boxsize_x;
			p.y  = ((double)rand()/(double)RAND_MAX-0.5)*r->boxsize_y;
			p.z  = 0.1*((double)rand()/(double)RAND_MAX-0.5)*r->boxsize_z;
			p.vx = 0;
			p.vy = -1.5*p.x*r->ri_sei.OMEGA;
			p.vz = 0;
			p.ax = 0; p.ay = 0; p.az = 0;
			p.m  = 0.01;
			p.r  = 0.05;
			particles_add(r, p);
		}
	}
	r->heartbeat = heartbeat;
	rebound_integrate(r,0);
}

void heartbeat(struct Rebound* const r){
	// Outputs a restartfile once per orbit.
	reb_output_timing(r,0);
	if (reb_output_check(r,1e-0*2.*M_PI)){
		reb_output_binary(r, "restart.bin");
		printf("\nSaved binary file.\n");
	}
}
