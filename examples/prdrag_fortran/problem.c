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
#include "tools.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "problem.h"
#include "integrator.h"


#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

/*extern double G;*/
/*extern int N;*/
/*extern int N_active;*/
void additional_forces();
const double k	 	= 0.01720209895;	// Gaussian constant 
double betaparticles = 0.01;    // beta parameter
                                // defined as the ratio of radiation pressure over gravity

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 		= 5;			// in days
	N_active	= 2;			// we treat pluto as a test particle. If all your particles have mass, remove this line.
	tmax		= 7.3e10;		// 200 Myr
	G		= k*k;			// These are the same units as used by the mercury6 code.
#ifdef OPENGL
	display_wire	= 1;			// Show orbits.
#endif // OPENGL
    problem_additional_forces = additional_forces;
	init_boxwidth(20); 			// Init box with width X astronomical units

	
	// Initial conditions for star
	struct particle star;
	star.x  = -4.06428567034226e-3;
    star.y  = -6.08813756435987e-3;
    star.z  = -1.66162304225834e-6; 
	star.vx = +6.69048890636161e-6;
    star.vy = -6.33922479583593e-6;
    star.vz = -3.13202145590767e-9;
	star.m  = 1.00000597682; 	// Sun + inner planets
	star.ax  = 0;
    star.ay  = 0;
    star.az  = 0; 
	particles_add(star);

	struct particle planet;
	planet.x  = +3.40546614227466e+0;
    planet.y  = +3.62978190075864e+0;
    planet.z  = +3.42386261766577e-2; 
	planet.vx = -5.59797969310664e-3;
    planet.vy = +5.51815399480116e-3;
    planet.vz = -2.66711392865591e-6;
	planet.m  = 1./1047.355;	// Jupiter
	planet.ax  = 0;
    planet.ay  = 0;
    planet.az  = 0; 
	particles_add(planet);

    // dust particles are initially on a circular orbit
    while(N<200){
        struct particle p;
        p.m  = 0;                   // massless
        double a = 1.;                  // a = 1 AU
        double v = sqrt(G*(star.m*(1.-betaparticles))/a);
        double phi = tools_uniform(0,2.*M_PI);      // random phase
        p.x  = a*sin(phi);  p.y  = a*cos(phi); p.z  = 0;
        p.vx = -v*cos(phi); p.vy = v*sin(phi); p.vz = 0;
        particles_add(p);
    }

	tools_move_to_center_of_momentum();
}

void problem_output(){
    if (output_check(10.*dt)){
        output_timing();
    }
}

void problem_finish(){
}

void problem_inloop(){
}


/*void for_sub(int nbod, struct particle *p);*/
void for_sub(int nbod);
struct particle p;

void additional_forces(){
    tidal_forces(&N, particles);
    force_radiation(&N, particles);
}
