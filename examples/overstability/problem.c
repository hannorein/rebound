/**
 * @file 	problem.c
 * @brief 	Example problem: overstability.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	A narrow box of Saturn's rings is simulated to
 * study the viscous overstability. Collisions are resolved using
 * the plane-sweep method.
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
#include "collisions.h"

extern double OMEGA;
extern double OMEGAZ;
extern double coefficient_of_restitution; 
extern double minimum_collision_velocity;

void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA 				= 1.;
	OMEGAZ 				= 3.6;
	dt				= 2e-3*2.*M_PI/OMEGA;
	double particle_r 		= 1;
	double tau			= 1.64;
	coefficient_of_restitution 	= 0.5;
	boxsize 			= 1;
	root_nx = 200; 	root_ny = 5; 	root_nz = 20;
	nghostx = 1; 	nghosty = 1; 	nghostz = 0;
	init_box();

	// Initial conditions
	double _N = tau * boxsize_x * boxsize_y/(M_PI*particle_r *particle_r);
	while (N<_N){
		struct particle p;
		p.x 	= ((double)rand()/(double)RAND_MAX-0.5)*boxsize_x;
		p.y 	= ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
		p.z 	= 10.0*((double)rand()/(double)RAND_MAX-0.5)*particle_r;
		p.vx 	= 0;
		p.vy 	= -1.5*p.x*OMEGA;
		p.vz 	= 0;
		p.ax 	= 0; p.ay 	= 0; p.az 	= 0;
		p.m 	= 1.;
		p.r 	= particle_r;
		particles_add(p);
	}
}

void problem_inloop(){
}

void problem_output(){
	if (output_check(2.*M_PI)){
		output_timing();
	}
}

void problem_finish(){
}
