/**
 * @file 	problem.c
 * @brief 	Example problem: bouncing balls at corner.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example tests collision detection methon near box boundaries.
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
	dt 		= 1e-3;
	boxsize 	= 3;
	coefficient_of_restitution = 1; // elastic collisions
	init_box();
	
	// Initial conditions
	int problem_id = 1;
	if (argc>1){			
		problem_id = atoi(argv[1]);
	}
	struct particle p;

	switch (problem_id){
		case 1: // Multiple instantaneous collisions across boundaries
			nghostx = 1; nghosty = 1; nghostz = 0;
			p.x  = 1; p.y  = 1; p.z  = 0;
			p.vx = 0; p.vy = 0; p.vz = 0;
			p.ax = 0; p.ay = 0; p.az = 0;
			p.m  = 1;
			p.r  = 0.1;
			particles_add(p);
			p.x  = -1; p.y  = -1; p.z  = 0;
			p.vx = 0;  p.vy = 0;  p.vz = 0;
			p.ax = 0;  p.ay = 0;  p.az = 0;
			p.m  = 1;
			p.r  = 0.1;
			particles_add(p);
			p.x  = 1; p.y  = -1; p.z  = 0;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 1;
			p.r  = 0.1;
			particles_add(p);
			p.x  = -1; p.y  = 1; p.z  = 0;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 1;
			p.r  = 0.1;
			particles_add(p);
			break;
		case 2: // Multiple instantaneous collisions with different sizes
			nghostx = 0; nghosty = 0; nghostz = 0;
			p.x  = 0; p.y  = 0; p.z  = 0;
			p.vx = 0; p.vy = 0; p.vz = 0;
			p.ax = 0; p.ay = 0; p.az = 0;
			p.m  = 1;
			p.r  = 0.5;
			particles_add(p);
			p.x  = 1; p.y  = 1; p.z  = 0;
			p.vx = 0; p.vy = 0; p.vz = 0;
			p.ax = 0; p.ay = 0; p.az = 0;
			p.m  = 0.008;
			p.r  = 0.1;
			particles_add(p);
			p.x  = -1; p.y  = -1; p.z  = 0;
			p.vx = 0;  p.vy = 0;  p.vz = 0;
			p.ax = 0;  p.ay = 0;  p.az = 0;
			p.m  = 0.008;
			p.r  = 0.1;
			particles_add(p);
			p.x  = 1; p.y  = -1; p.z  = 0;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 0.008;
			p.r  = 0.3;
			particles_add(p);
			p.x  = -1; p.y  = 1; p.z  = 0;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 0.008;
			p.r  = 0.2;
			particles_add(p);
			p.x  = 0;  p.y  = 0; p.z  = 1.3;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 0.008;
			p.r  = 0.1;
			particles_add(p);
			p.x  = 0;  p.y  = 0; p.z  =-1.3;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 0.008;
			p.r  = 0.05;
			particles_add(p);
			break;
	}
}

void problem_inloop(){
}

void problem_output(){
}

void problem_finish(){
}
