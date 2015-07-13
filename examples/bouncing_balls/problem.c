/**
 * @file 	problem.c
 * @brief 	Example problem: bouncing balls.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example is a simple test of collision detection 
 * methods. To change the collision detection algorithm, you can replace
 * the module collisions_direct.c to either collisions_tree.c or 
 * collisions_sweep.c in the Makefile. 
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
#include "particle.h"
#include "rebound.h"

extern double coefficient_of_restitution; 

int main(int argc, char* argv[]){
	struct Rebound* r = rebound_init();
	// Setup constants
	r->integrator	= LEAPFROG;
	r->gravity	= RB_GT_BASIC;
	r->collision	= RB_CT_DIRECT;
	r->dt = 1e-2;
	rebound_configure_box(r, 3.0,1,1,1);
	coefficient_of_restitution = 1; // elastic collisions
	
	// Initial conditions
	{
		struct Particle p;
		p.x  = 1; p.y  = 1; p.z  = 1;
		p.vx = 0; p.vy = 0; p.vz = 0;
		p.ax = 0; p.ay = 0; p.az = 0;
		p.m  = 1;
		p.r  = 0.1;
		particles_add(r, p);
	}
	{
		struct Particle p;
		p.x  = -1; p.y  = -1; p.z  = -1;
		p.vx =  0; p.vy =  0; p.vz =  0;
		p.ax =  0; p.ay =  0; p.az =  0;
		p.m  = 1;
		p.r  = 0.1;
		particles_add(r, p);
	}

	rebound_integrate(r, 0);
}

