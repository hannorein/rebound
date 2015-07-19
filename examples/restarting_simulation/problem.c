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

void heartbeat(struct reb_simulation* const r);

int main(int argc, char* argv[]){
	{
		printf("Running simulation until t=10.\n");
		struct reb_simulation* r = reb_create_simulation();
		r->integrator	= REB_INTEGRATOR_SEI;
		r->collision	= REB_COLLISION_DIRECT;
		r->ri_sei.OMEGA	= 1.;	
		r->dt 		= 1e-4*2.*M_PI; 
		r->nghostx = 1; r->nghosty = 1; r->nghostz = 0;
		reb_configure_box(r,1.,1,1,1);

		// Check if simulation is restarted
		// Initial conditions
		while (r->N<50){
			struct reb_particle p;
			p.x  = ((double)rand()/(double)RAND_MAX-0.5)*r->boxsize.x;
			p.y  = ((double)rand()/(double)RAND_MAX-0.5)*r->boxsize.y;
			p.z  = 0.1*((double)rand()/(double)RAND_MAX-0.5)*r->boxsize.z;
			p.vx = 0;
			p.vy = -1.5*p.x*r->ri_sei.OMEGA;
			p.vz = 0;
			p.ax = 0; p.ay = 0; p.az = 0;
			p.m  = 0.01;
			p.r  = 0.05;
			reb_add(r, p);
		}
		reb_integrate(r,10.);
		printf("Saving simulation to binary file and freeing up memory.\n");
		reb_output_binary(r, "restart.bin");
		reb_free_simulation(r);
		r = NULL;
	}
	{
		printf("Creating simulation from binary file and integrating until t=20.\n");
		struct reb_simulation* r = reb_create_simulation_from_binary("restart.bin");
		reb_integrate(r,20.);
		printf("Done.\n");
	}
	
}

