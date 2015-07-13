/**
 * @file 	problem.c
 * @brief 	Example problem: bouncing string.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example tests collision detection methods.
 * The example uses a non-square, rectangular box. 10 particles are placed
 * along a line. All except one of the particles are at rest 
 * initially.
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
#include "particle.h"

extern double coefficient_of_restitution; 
int main(int argc, char* argv[]){
	struct Rebound* const r = rebound_init();
	// Setup constants
	r->dt 				= 1e-3;
	rebound_configure_box(r,10.,3,1,1);  // boxsize 10., three root boxes in x direction, one in y and z
	coefficient_of_restitution 	= 1; // elastic collisions
	r->integrator			= LEAPFROG;
	r->boundary			= RB_BT_PERIODIC;
	r->nghostx = 1; 
	r->nghosty = 1; 
	r->nghostz = 0;

	// Initial conditions
	for(int i=0;i<10;i++){
		struct Particle p;
		p.x  = -r->boxsize_x/2.+r->boxsize_x*(double)i/10.;
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
		particles_add(r, p);
	}

	// Give one particle a kick
	r->particles[0].vx = 20;

	rebound_integrate(r,0);
}
