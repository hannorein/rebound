/**
 * @file 	problem.c
 * @brief 	Example problem: spreading ring.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	A narrow ring is spreading due to collisions.
 * We use the Wisdom Holman integrator and a plane-sweep algorithm 
 * in the phi direction.
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
#include "collisions.h"
#include "tools.h"

extern double opening_angle2;
extern int Nmax;

void problem_init(int argc, char* argv[]){
	// Setup constants
	G 		= 1;		
	N_active	= 1;
	N_collisions	= 1; 	// Don't detect collisions for the star.
	softening 	= 0.01;		
	dt 		= 1e-3;
	boxsize 	= 2.4;
	root_nx = 1; root_ny = 1; root_nz = 1;
	nghostx = 0; nghosty = 0; nghostz = 0; 		
	init_box();

	// Setup particles
	int _N = 1000;
	// Initial conditions
	struct particle star;
	star.x 		= 0; star.y 	= 0; star.z	= 0;
	star.vx 	= 0; star.vy 	= 0; star.vz 	= 0;
	star.ax 	= 0; star.ay 	= 0; star.az 	= 0;
	star.m 		= 1;
	star.r		= 0.01;
	particles_add(star);
	while(N<_N){
		struct particle pt;
		double a	= tools_powerlaw(boxsize/2.9,boxsize/3.1,.5);
		double phi 	= tools_uniform(0,2.*M_PI);
		pt.x 		= a*cos(phi);
		pt.y 		= a*sin(phi);
		pt.z 		= a*tools_normal(0.0001);
		double vkep 	= sqrt(G*star.m/a);
		pt.vx 		=  vkep * sin(phi);
		pt.vy 		= -vkep * cos(phi);
		pt.vz 		= 0;
		pt.ax 		= 0;
		pt.ay 		= 0;
		pt.az 		= 0;
		pt.m 		= 0.0001;
		pt.r 		= .3/sqrt((double)_N);
		particles_add(pt);
	}
}

void problem_inloop(){
	for (int i=N_active;i<N;i++){
		struct particle p = particles[i];
		double r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
		// Remove particles falling into the star.
		if (r<0.03){
			particles[i] = particles[N-1];
			i--;
			N--;
		}
	}
}

void problem_output(){
	if (output_check(10.0*dt)) output_timing();
}

void problem_finish(){
}
