/**
 * @file 	problem.c
 * @brief 	Example problem: restarting simulation.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example demonstrates how to restart a simulation
 * using a binary file.
 * 
 * First, run the program with './nbody'.
 * Random initial conditions are created and
 * a restart file is written once per orbit.
 * Then, to restart the simulation, run the 
 * program with './nbody --restart restart.bin'.
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
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "input.h"


extern double OMEGA;
extern double coefficient_of_restitution;
extern double minimum_collision_velocity;

void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA 		= 1.;
	dt 		= 1e-4*2.*M_PI; 
	boxsize 	= 1;
	coefficient_of_restitution = 0.5;
	minimum_collision_velocity = 0.01;
	nghostx = 1; nghosty = 1; nghostz = 0;
	init_box();

	// Check if simulation is restarted
	if (input_check_restart(argc,argv)!=1){
		// Initial conditions
		while (N<50){
			struct particle p;
			p.x  = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_x;
			p.y  = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
			p.z  = 0.1*((double)rand()/(double)RAND_MAX-0.5)*boxsize_z;
			p.vx = 0;
			p.vy = -1.5*p.x*OMEGA;
			p.vz = 0;
			p.ax = 0; p.ay = 0; p.az = 0;
			p.m  = 0.01;
			p.r  = 0.05;
			particles_add(p);
		}
	}
}

void problem_inloop(){
}

void problem_output(){
	// Outputs a restartfile once per orbit.
	output_timing();
	if (output_check(1e-0*2.*M_PI&&t!=0)){
		output_binary("restart.bin");
		printf("\nSaved binary file. Restart simulation with './nbody --restart restart.bin'.\n");
	}
}

void problem_finish(){
}
