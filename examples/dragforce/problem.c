/**
 * @file 	problem.c
 * @brief 	Example problem: Velocity dependent drag force
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This is a very simple example on how to implement a velocity 
 * dependent drag force. The example uses the IAS15 integrator, which 
 * is ideally suited to handle non-conservative forces.
 * No gravitational forces or collisions are present.
 * 
 * @section 	LICENSE
 * Copyright (c) 2013 Hanno Rein, Dave Spiegel
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
#include "tools.h"
#include "output.h"
#include "particle.h"
#include "problem.h"
#include "integrator.h"

void additional_forces();

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 			= 1e-4;		// initial timestep.
	boxsize 		= 10;	
	tmax			= 40;

	// Setup callback function for velocity dependent forces.
	problem_additional_forces 	= additional_forces;
	
	// Initialize simulation and particles
	init_box();
	
	
	struct particle p; 
	p.m  = 0;	// massless
	p.x = 1; 	p.y = 0; 	p.z = 0;
	p.vx = -1; 	p.vy = 0; 	p.vz = 0;
	p.ax = 0; 	p.ay = 0; 	p.az = 0;
	particles_add(p); 

	// Delete previous output
	system("rm -v r.txt");	
}

void additional_forces(){
	// Simplest velocity dependent drag force.
	double dragcoefficient = 1;
	for (int i=0;i<N;i++){
		particles[i].ax = -dragcoefficient*particles[i].vx;
		particles[i].ay = -dragcoefficient*particles[i].vy;
		particles[i].az = -dragcoefficient*particles[i].vz;
	}
}

void problem_output(){
	// Output some information to the screen every 100th timestep
	if(output_check(100.*dt)){
		output_timing();
	}
	// Output the particle position to a file every timestep.
	FILE* f = fopen("r.txt","a");
	fprintf(f,"%e\t%e\t%e\n",t,particles[0].x, particles[1].vx);
	fclose(f);
}

void problem_finish(){
}
