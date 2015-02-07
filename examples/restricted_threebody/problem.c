/**
 * @file 	problem.c
 * @brief 	Example problem: Restricted three body problem.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example simulates a disk of test particles around 
 * a central object, being perturbed by a planet. 
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

void problem_init(int argc, char* argv[]){
	// Setup constants
	boxsize 		= 8; 
	softening		= 1e-6;
	dt 			= 1.0e-2*2.*M_PI;
	N_active 		= 2; 	// Only the star and the planet have non-zero mass
	init_box();
	
	// Initial conditions for star
	struct particle star;
	star.x  = 0; 	star.y  = 0; 	star.z  = 0; 
	star.vx = 0; 	star.vy = 0; 	star.vz = 0;
	star.m  = 1;
	particles_add(star);

	// Initial conditions for planet
	double planet_e = 0.;
	struct particle planet;
	planet.x  = 1.-planet_e; 	planet.y  = 0; 				planet.z  = 0; 
	planet.vx = 0; 			planet.vy = sqrt(2./(1.-planet_e)-1.); 	planet.vz = 0;
	planet.m  = 1e-2;
	particles_add(planet);
	
	while(N<10000){
		double x 	= ((double)rand()/(double)RAND_MAX-0.5)*boxsize*0.9;
		double y 	= ((double)rand()/(double)RAND_MAX-0.5)*boxsize*0.9;
		double a 	= sqrt(x*x+y*y);
		double phi 	= atan2(y,x);
		if (a<.1) continue;
		if (a>boxsize_x/2.*0.9) continue;

		double vkep = sqrt(G*star.m/a);
		struct particle testparticle;
		testparticle.x  = x;
		testparticle.y  = y; 
		testparticle.z  = 1.0e-2*x*((double)rand()/(double)RAND_MAX-0.5);
		testparticle.vx = -vkep*sin(phi);
		testparticle.vy = vkep*cos(phi);
		testparticle.vz = 0;
		testparticle.ax = 0; 
		testparticle.ay = 0; 
		testparticle.az = 0;
		testparticle.m  = 0;
		particles_add(testparticle);
	}
}

void problem_output(){
	output_timing();
	if (output_check(2.*M_PI)){
		output_orbits("orbit.txt");
	}
}

void problem_finish(){
}

