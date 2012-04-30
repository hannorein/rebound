/**
 * @file 	problem.c
 * @brief 	Example problem: circular orbit.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example uses the Wisdom Holman integrator
 * to integrate particles on a circular orbit in a fixed 
 * potential.
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

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 		= M_PI*1e-3*1.1234125235345;
	boxsize 	= 2.8;
	N_active      = 2;
	init_box();

	// Initial conditions
	// The WH integrator assumes a heliocentric coordinate system.
	// Therefore the star has to be at the origin.
	struct particle star;
	star.x  = 0; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = 0; star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = 1;
	particles_add(star);

	// Planet 1 at 1 AU
	struct particle planet;
	planet.m  = 1.e-3;
	planet.x  = 1.; planet.y  = 0.; planet.z  = 0.;
	planet.ax = 0; planet.ay = 0; planet.az = 0;
	planet.vx = 0;
	planet.vy = sqrt(G*(star.m+planet.m)/planet.x);
	planet.vz = 0;
	particles_add(planet);

	// Ring particles
	double inner_R =  1./40.;
	double outer_R =  2./40.;
	int _N = 5000;
	for (int i = 0; i < _N; i++) {
		double dr = inner_R + (double)i/(double)(_N-1) * (outer_R-inner_R);
		double dvr = sqrt(planet.m / dr);
		struct particle ringparticle = planet;
		ringparticle.x  += dr;
		ringparticle.vy += dvr;
		ringparticle.m  = 0;
		particles_add(ringparticle);
	}
}

void problem_inloop(){
  if(output_check(1000.*dt)){
    output_timing();
  }
  if(output_check(124.)){
    //output_append_ascii("position.txt");
    output_append_orbits("orbits.txt");
  }
}

void problem_output(){
}

void problem_finish(){
}
