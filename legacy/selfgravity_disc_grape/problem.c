/**
 * @file 	problem.c
 * @brief 	Example problem: self-gravity disc.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	A self-gravitating disc is integrated using
 * the leap frog integrator. This example is using the GRAPE
 * module to calculate the self-gravity. You need to have a physical
 * GRAPE card in your computer to run this example.
 * Collisions are not resolved.
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

extern int Nmax;

int main(int argc, char* argv[]){
	// Setup constants
	G 		= 1;		
	integrator	= LEAPFROG;
	softening 	= 0.01;		
	dt 		= 3e-3;
	boxsize 	= 1.2;
	N_root_x = 1; N_root_y = 1; N_root_z = 1;
	N_ghost_x = 0; N_ghost_y = 0; N_ghost_z = 0; 		
	init_box();

	// Setup particles
	double disc_mass = 2e-1;
	int _N = 10000;
	// Initial conditions
	struct reb_particle star;
	star.x 		= 0; star.y 	= 0; star.z	= 0;
	star.vx 	= 0; star.vy 	= 0; star.vz 	= 0;
	star.ax 	= 0; star.ay 	= 0; star.az 	= 0;
	star.m 		= 1;
	reb_simulation_add(r, star);
	while(N<_N){
		struct reb_particle pt;
		double a	= reb_random_powerlaw(boxsize/10.,boxsize/2./1.2,-1.5);
		double phi 	= reb_random_uniform(0,2.*M_PI);
		pt.x 		= a*cos(phi);
		pt.y 		= a*sin(phi);
		pt.z 		= a*reb_random_normal(0.001);
		double mu 	= star.m + disc_mass * (pow(a,-3./2.)-pow(boxsize/10.,-3./2.))/(pow(boxsize/2./1.2,-3./2.)-pow(boxsize/10.,-3./2.));
		double vkep 	= sqrt(G*mu/a);
		pt.vx 		=  vkep * sin(phi);
		pt.vy 		= -vkep * cos(phi);
		pt.vz 		= 0;
		pt.ax 		= 0;
		pt.ay 		= 0;
		pt.az 		= 0;
		pt.m 		= disc_mass/(double)_N;
		reb_simulation_add(r, pt);
	}
}

void heartbeat(struct reb_simulation* r){
	if (reb_simulation_output_check(10.0*dt)) reb_simulation_output_timing();
	if (reb_simulation_output_check(1.)) reb_simulation_output_ascii("ascii.txt");
}

void problem_finish(){
}
