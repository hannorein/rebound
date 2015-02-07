/**
 * @file 	problem.c
 * @brief 	Example problem: Close Encounter.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This example integrates a densly packed planetary system 
 * which becomes unstable on a timescale of only a few orbits. The IAS15 
 * integrator with adaptive timestepping is used. This integrator 
 * automatically decreases the timestep whenever a close 
 * enocunter happens. IAS15 is very high order and ideally suited for the 
 * detection of these kind of encounters.
 * 
 * @section 	LICENSE
 * Copyright (c) 2013 Hanno Rein, Dave Spiegel
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
#include "tools.h"
#include "output.h"
#include "particle.h"
#include "integrator.h"

#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

void problem_init(int argc, char* argv[]){
	dt = 0.01*2.*M_PI;						// initial timestep
	// integrator_epsilon = 1e-2;					// accuracy parameter, default is 1e-2 and should work in most cases.

#ifdef OPENGL
	display_wire	= 1;						// show instantaneous orbits
#endif // OPENGL
	init_boxwidth(10); 					

	struct particle star;
	star.m = 1;
	star.x = 0; 	star.y = 0; 	star.z = 0;
	star.vx = 0; 	star.vy = 0; 	star.vz = 0;
	particles_add(star);
	
	// Add planets
	int N_planets = 7;
	for (int i=0;i<N_planets;i++){
		double a = 1.+(double)i/(double)(N_planets-1);		// semi major axis
		double v = sqrt(1./a); 					// velocity (circular orbit)
		struct particle planet;
		planet.m = 1e-4; 
		planet.x = a; 	planet.y = 0; 	planet.z = 0;
		planet.vx = 0; 	planet.vy = v; 	planet.vz = 0;
		particles_add(planet); 
	}
	tools_move_to_center_of_momentum();				// This makes sure the planetary systems stays within the computational domain and doesn't drift.
}

void problem_output(){
	if (output_check(10.*2.*M_PI)){  
		output_timing();
	}
}

void problem_finish(){
}
