/**
 * @file 	problem.c
 * @brief 	Example problem: Close Encounter with output.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This example integrates a densly packed planetary system 
 * which becomes unstable on a timescale of only a few orbits. 
 * The example is identical to the `close_encounter` sample, except that 
 * the collisions are recorded and written to a file. What kind of collisions
 * are recorded can be easily modified. It is also possible to implement some
 * additional physics whenever a collision has been detection (e.g. fragmentation).
 * The collision search is by default a direct search, i.e. O(N^2) but can be
 * changed to a tree by using the `collisions_tree.c` module.
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
#include "collision_resolve.h"

#ifdef OPENGL
extern int display_wire;
extern int display_spheres;
#endif // OPENGL

// Define our own collision resolve function, which will only record collisions but not change any of the particles.		
void collision_record_only(struct collision c){
	double delta_t = 2.*M_PI; 					// 1 year	
	
	// only record a maximum of one collision per year per particle
	if ( particles[c.p1].lastcollision+delta_t < t  &&  particles[c.p2].lastcollision+delta_t < t ){
		particles[c.p1].lastcollision = t; 
		particles[c.p2].lastcollision = t;
		FILE* of = fopen("collisions.txt","a+");		// open file for collision output
		fprintf(of, "%e\t",t);					// time
		fprintf(of, "%e\t",(particles[c.p1].x+particles[c.p2].x)/2.);	// x position
		fprintf(of, "%e\t",(particles[c.p1].y+particles[c.p2].y)/2.);	// y position
		fprintf(of, "\n");
		fclose(of);						// close file
	}
}


void problem_init(int argc, char* argv[]){
	dt = 0.1*2.*M_PI;						// initial timestep
	// integrator_epsilon = 1e-2;					// accuracy parameter, default is 1e-2 and should work in most cases.

#ifdef OPENGL
	display_wire	= 1;						// show instantaneous orbits
	display_spheres = 0;						// don't show spheres
#endif // OPENGL
	init_boxwidth(10); 					

	collision_resolve = collision_record_only;			// Set function pointer for collision recording.

	struct particle star;
	star.m = 1;
	star.r = 0;							// Star is pointmass
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
		double rhill = a * pow(planet.m/(3.*star.m),1./3.);	// Hill radius
		planet.r = rhill;					// Set planet radius to hill radius 
									// A collision is recorded when planets get within their hill radius
									// The hill radius of the particles might change, so it should be recalculated after a while
		planet.lastcollision = 0; 
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
