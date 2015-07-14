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

double coefficient_of_restitution(const struct reb_context* const r, double vel){
	return 1.;	
} 
int main(int argc, char* argv[]){
	struct reb_context* const r = reb_init();
	// Setup constants
	r->dt 				= 1e-3;
	reb_configure_box(r,10.,3,1,1);  // boxsize 10., three root boxes in x direction, one in y and z
	r->coefficient_of_restitution 	= coefficient_of_restitution; // elastic collisions
	r->integrator			= RB_IT_LEAPFROG;
	r->boundary			= RB_BT_PERIODIC;
	r->collision			= RB_CT_DIRECT;
	r->gravity			= RB_GT_NONE;
	r->nghostx = 1; 
	r->nghosty = 1; 
	r->nghostz = 0;

	// Initial conditions
	for(int i=0;i<10;i++){
		struct reb_particle p;
		p.x  = -r->boxsize.x/2.+r->boxsize.x*(double)i/10.;
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
		reb_add(r, p);
	}

	// Give one particle a kick
	r->particles[0].vx = 20;

	reb_integrate(r,0);
}
