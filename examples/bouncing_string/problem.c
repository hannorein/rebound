/**
 * @file 	problem.c
 * @brief 	Example problem: bouncing string.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example tests collision detection methods.
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

extern double coefficient_of_restitution; 
void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 				= 1e-3;
	tmax 				= 10000;
	boxsize 			= 10;
	coefficient_of_restitution 	= 1; // elastic collisions
	root_nx = 3; root_ny = 1; root_nz = 1;
	nghostx = 1; nghosty = 1; nghostz = 0;
	init_box();

	// Initial conditions
	while(N<10){
		struct particle p;
		p.x  = -boxsize_x/2.+boxsize_x*(double)N/10.;
		p.y  = 0;
		p.z  = 0;
		p.vx = 0;
		p.vy = 0;
		p.vz = 0;
		p.ax = 0;
		p.ay = 0;
		p.az = 0;
		p.m  = 1;
		p.r  = 1;
		particles_add(p);
	}
	particles[0].vx = 20;
}

void problem_inloop(){
}

void problem_output(){
}

void problem_finish(){
}
