/**
 * @file 	problem.c
 * @brief 	Example problem: Kozai.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example uses the IAS15 integrator to simulate
 * a very eccentric planetary orbit. The integrator
 * automatically adjusts the timestep so that the pericentre passages
 * resovled with high accuracy.
 *
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
	// Setup constants
	G			= 1;		// Gravitational constant
#ifdef OPENGL
	display_wire		= 1; 		// show istantaneous orbits.
#endif // OPENGL


	double e_testparticle 	= 1.-1e-7;	
	double mass_scale	= 1.;		// Some integrators have problems when changing the mass scale, IAS15 does not. 
	double size_scale	= 1;		// Some integrators have problems when changing the size scale, IAS15 does not.


	boxsize 		= 25.*size_scale;
	init_box();
	
	struct particle star; 
	star.m  = mass_scale;
	star.x  = 0; star.y  = 0; star.z  = 0; 
	star.vx = 0; star.vy = 0; star.vz = 0;
	particles_add(star); 
	
	struct particle planet; 
	planet.m  = 0;
	planet.x  = size_scale*(1.-e_testparticle); planet.y  = 0; planet.z  = 0; 
	planet.vx = 0; planet.vy = sqrt((1.+e_testparticle)/(1.-e_testparticle)*mass_scale/size_scale); planet.vz = 0;
	particles_add(planet); 
	
	tools_move_to_center_of_momentum();
	
	// initial timestep
	dt 			= 1e-13*sqrt(size_scale*size_scale*size_scale/mass_scale); 
	tmax			= 1e2*2.*M_PI*sqrt(size_scale*size_scale*size_scale/mass_scale);
	
}

void problem_output(){
	if(output_check(tmax/10000.)){		// outputs to the screen
		output_timing();
	}
	// Output the time and the current timestep. Plot it to see how IAS15 automatically reduces the timestep at pericentre. 
	FILE* of = fopen("timestep.txt","a"); 
	fprintf(of,"%e\t%e\t\n",t/tmax,dt/tmax);
	fclose(of);
}

void problem_finish(){
}
