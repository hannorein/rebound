/**
 * @file 	problem.c
 * @brief 	Example problem: Kozai.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
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
#include "output.h"
#include "particle.h"
#include "boundaries.h"

#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 		= M_PI*1e-2*1.1234125235345; // The random number ensures that the time step is not a multiple of any other frequency.
	boxsize 	= 25;
	tmax		= 1.6e4;
#ifdef OPENGL
	display_wire	= 1; // Show istantaneous orbits..
#endif // OPENGL
	init_box();

	// Initial conditions
	struct particle p; 
	
	// Star
	p.m  = 1;
	// The WH integrator assumes a heliocentric coordinate system. 
	// Therefore the star has to be at the origin. 
	p.x  = 0; p.y  = 0; p.z  = 0; 
	p.vx = 0; p.vy = 0; p.vz = 0;
	p.ax = 0; p.ay = 0; p.az = 0;
	particles_add(p); 
	
	// Test particle
	// Actually this is treated as a 'massive particle with zero mass'.
	// This ensures correct ordering of Jacobi coordinates in the WH integrator. 
	p.m  = 0;
	double e_testparticle = 0;
	p.x  = 1.-e_testparticle; p.y  = 0; p.z  = 0; 
	p.vx = 0; p.vy = sqrt((1.+e_testparticle)/(1.-e_testparticle)); p.vz = 0;
	p.ax = 0; p.ay = 0; p.az = 0;
	particles_add(p); 
	
	// Perturber
	p.x  = 10; p.y  = 0; p.z  = 0; 
	double inc_perturber = 89.9;
	p.vx = 0; 
	p.m  = 1;
	p.vy = cos(inc_perturber/180.*M_PI)*sqrt((1.+p.m)/p.x); 
	p.vz = sin(inc_perturber/180.*M_PI)*sqrt((1.+p.m)/p.x); 
	p.ax = 0; p.ay = 0; p.az = 0;
	particles_add(p); 
	
	system("rm -v orbits.txt");	
}

void problem_inloop(){
	if(output_check(4000.*dt)){
		output_timing();
	}
	if(output_check(12.)){
		output_append_orbits("orbits.txt");
	}
}

void problem_output(){
}

void problem_finish(){
}
