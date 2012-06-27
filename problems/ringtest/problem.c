/**
 * @file 	problem.c
 * @brief 	Example problem: circular orbit.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example uses the Wisdom Holman integrator
 * to integrate particles on a circular orbit in a fixed 
 * potential.
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
#include <string.h>
#include "main.h"
#include "problem.h"
#include "input.h"
#include "tools.h"
#include "output.h"
#include "particle.h"
#include "boundaries.h"

// Star and planet (note, those wont be updated after they have been inserted)
#ifdef INTEGRATOR_RADAU15
extern double 	integrator_accuracy;	// Desired accuracy. Play with this, make sure you get a converged results.
extern double 	integrator_min_dt;	// Timestep floor.
extern int 	integrator_adaptive_timestep;
#endif // INTEGRATOR_RADAU15

double massfac = 1;
double distancefac = 1;

void problem_init(int argc, char* argv[]){
	// Setup constants
	G				= 1;
	dt 				= 1e-3;
#ifdef INTEGRATOR_RADAU15
	integrator_accuracy 		= 1e1;
	integrator_adaptive_timestep	= 1;
#endif // INTEGRATOR_RADAU15
	boxsize 			= 10.*distancefac;
	N_active      			= 2;
	init_box();

	// Initial conditions
	struct particle star;
	star.x  = 0.5*distancefac; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = .5*sqrt(massfac)/sqrt(distancefac); star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = .5*massfac;
	particles_add(star);
	
	star.x  = -0.5*distancefac; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = -.5*sqrt(massfac)/sqrt(distancefac); star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = .5*massfac;
	particles_add(star);

//	tools_move_to_center_of_momentum();
}
void output_dt(){
	FILE* of = fopen("dt.txt","a+"); 
	double period = 2.*M_PI*sqrt(distancefac*distancefac*distancefac/massfac);
	fprintf(of,"%e\t%e\n",t/period,dt/period);
	fclose(of);
}
		
void problem_output(){
	// Output stuff
	if(output_check(2.*M_PI)){
		output_timing();
	}
	output_dt();
}
void problem_inloop(){
}

void problem_finish(){
}
