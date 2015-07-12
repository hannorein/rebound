/**
 * @file 	problem.c
 * @brief 	Example problem: bouncing balls at corner.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example tests collision detection methods across box boundaries.
 * There are four particles, one in each corner. To see the ghost boxes in OpenGL
 * press `g` while the simulation is running.
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
	struct Rebound* r = rebound_init();
	// Setup constants
	r->integrator	= LEAPFROG;
	r->dt 		= 1e-3;
	r->boundary	= RB_BT_PERIODIC;
	rebound_configure_box(r,3.,1,1,1);
	coefficient_of_restitution = 1; // elastic collisions
	
	// Initial conditions
	int problem_id = 1;
	if (argc>1){			
		problem_id = atoi(argv[1]);
	}
	struct Particle p;

	switch (problem_id){
		case 1: // Multiple instantaneous collisions across boundaries
			r->nghostx = 1; r->nghosty = 1; r->nghostz = 0;
			p.x  = 1; p.y  = 1; p.z  = 0;
			p.vx = 0; p.vy = 0; p.vz = 0;
			p.ax = 0; p.ay = 0; p.az = 0;
			p.m  = 1;
			p.r  = 0.1;
			particles_add(r,p);
			p.x  = -1; p.y  = -1; p.z  = 0;
			p.vx = 0;  p.vy = 0;  p.vz = 0;
			p.ax = 0;  p.ay = 0;  p.az = 0;
			p.m  = 1;
			p.r  = 0.1;
			particles_add(r,p);
			p.x  = 1; p.y  = -1; p.z  = 0;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 1;
			p.r  = 0.1;
			particles_add(r,p);
			p.x  = -1; p.y  = 1; p.z  = 0;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 1;
			p.r  = 0.1;
			particles_add(r,p);
			break;
		case 2: // Multiple instantaneous collisions with different sizes
			r->nghostx = 0; r->nghosty = 0; r->nghostz = 0;
			p.x  = 0; p.y  = 0; p.z  = 0;
			p.vx = 0; p.vy = 0; p.vz = 0;
			p.ax = 0; p.ay = 0; p.az = 0;
			p.m  = 1;
			p.r  = 0.5;
			particles_add(r,p);
			p.x  = 1; p.y  = 1; p.z  = 0;
			p.vx = 0; p.vy = 0; p.vz = 0;
			p.ax = 0; p.ay = 0; p.az = 0;
			p.m  = 0.008;
			p.r  = 0.1;
			particles_add(r,p);
			p.x  = -1; p.y  = -1; p.z  = 0;
			p.vx = 0;  p.vy = 0;  p.vz = 0;
			p.ax = 0;  p.ay = 0;  p.az = 0;
			p.m  = 0.008;
			p.r  = 0.1;
			particles_add(r,p);
			p.x  = 1; p.y  = -1; p.z  = 0;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 0.008;
			p.r  = 0.3;
			particles_add(r,p);
			p.x  = -1; p.y  = 1; p.z  = 0;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 0.008;
			p.r  = 0.2;
			particles_add(r,p);
			p.x  = 0;  p.y  = 0; p.z  = 1.3;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 0.008;
			p.r  = 0.1;
			particles_add(r,p);
			p.x  = 0;  p.y  = 0; p.z  =-1.3;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 0.008;
			p.r  = 0.05;
			particles_add(r,p);
			break;
	}

	rebound_integrate(r,0);
}

