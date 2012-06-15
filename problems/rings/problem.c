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
#include <string.h>
#include "main.h"
#include "problem.h"
#include "input.h"
#include "tools.h"
#include "output.h"
#include "particle.h"
#include "boundaries.h"

// Star and planet (note, those wont be updated after they have been inserted)
struct particle star;
struct particle planet;
extern double 	integrator_accuracy;	// Desired accuracy. Play with this, make sure you get a converged results.
extern double 	integrator_min_dt;	// Timestep floor.
void radiation_forces();

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 				= 1e-1;
	// integrator_min_dt		= 1e-2;		// Set a timestep floor.
	integrator_accuracy 		= 1e-4;
	boxsize 			= 2.8;
	N_active      			= 2;
	problem_additional_forces 	= radiation_forces;
	init_box();

	// Initial conditions
	star.x  = 0; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = 0; star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = 1;
	particles_add(star);

	// Planet 1 at 1 AU
	planet.m  = 1.e-3;
	planet.x  = 1.; planet.y  = 0.; planet.z  = 0.;
	planet.ax = 0; planet.ay = 0; planet.az = 0;
	planet.vx = 0;
	planet.vy = sqrt(G*(star.m+planet.m)/planet.x);
	planet.vz = 0;
	particles_add(planet);
	output_double("planet mass", planet.m);
	output_double("planet a", planet.x);

	// Ring particles
	double planet_hill = planet.x*pow(planet.m/3./star.m,1./3.);
	output_double("planet hill", planet_hill);
	double r_inner =  planet_hill*input_get_double(argc, argv, "r_inner", 0.1);
	double r_outer =  planet_hill*input_get_double(argc, argv, "r_outer", 0.4);
	output_double("ring inner", r_inner);
	output_double("ring outer", r_outer);
	int _N = 1000;
	for (int i = 0; i < _N; i++) {
		struct particle ringparticle = planet;
		ringparticle.m  = 0;
		double r 	= tools_uniform(r_inner,r_outer);
		double v	= sqrt(G*planet.m / r);
		double phi	= tools_uniform(0,2.*M_PI);
		ringparticle.x  +=  r*cos(phi);
		ringparticle.y  +=  r*sin(phi);
		ringparticle.vx += -v*sin(phi);
		ringparticle.vy +=  v*cos(phi);
		particles_add(ringparticle);
	}
	system("cat config.log");
	tools_move_to_center_of_momentum();
}

void radiation_forces(){
	// Star 
	const struct particle star = particles[0];				// cache
#pragma omp parallel for
	for (int i=2;i<N;i++){
		const struct particle p = particles[i]; 			// cache
		const double prx  = p.x-star.x;
		const double pry  = p.y-star.y;
		const double prz  = p.z-star.z;
		const double prvx = p.vx-star.vx;
		const double prvy = p.vy-star.vy;
		const double prvz = p.vz-star.vz;

		const double beta	= 0.012;
		const double pr 	= sqrt(prx*prx + pry*pry + prz*prz); 	// distance relative to star
		const double c 		= 10064.915; 				// speed of light.
		const double rdot 	= (prvx*prx + prvy*pry + prvz*prz)/pr; 	// radial velocity relative to  star
		const double F_r 	= beta*G*star.m/(pr*pr);

		particles[i].ax += F_r*((1.-rdot/c)*prx/pr - prvx/c);
		particles[i].ay += F_r*((1.-rdot/c)*pry/pr - prvy/c);
		particles[i].az += F_r*((1.-rdot/c)*prz/pr - prvz/c);
	}
}


void problem_output(){
	if(output_check(2.*M_PI)){
		output_timing();
	}
	if(output_check(124.)){
		//output_append_ascii("position.txt");
		//output_append_orbits("orbits.txt");
	}
}
void problem_inloop(){
}

void problem_finish(){
}
