/**
 * @file 	problem.c
 * @brief 	Example problem: forced migration of GJ876.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example applies dissipative forces to two
 * bodies orbiting a central object. The forces are specified
 * in terms of damping timescales for the semi-major axis and
 * eccentricity. This mimics planetary micration in a proto-
 * stellar disc. The example reproduces the study of Lee & 
 * Peale (2002) on the formation of the planetary system 
 * GJ876. For a comparison, see figure 4 in their paper.
 *
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
#include "output.h"
#include "particle.h"
#include "boundaries.h"

double* tau_a; 	/**< Migration timescale in years for all particles */
double* tau_e; 	/**< Eccentricity damping timescale in years for all particles */

#ifdef OPENGL
extern int display_wire;
#endif 	// OPENGL

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 		= 1e-3*2.*M_PI;
	boxsize 	= 3;
#ifdef OPENGL
	display_wire 	= 1;
#endif 	// OPENGL
	init_box();

	// Initial conditions
	// Parameters are those of Lee & Peale 2002, Figure 4. 
	struct particle star;
	star.x  = 0; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = 0; star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = 0.32;		// This is a sub-solar mass star
	particles_add(star); 
	
	struct particle p1;	// Planet 1
	p1.x 	= 0.5;	p1.y = 0;	p1.z = 0;
	p1.ax 	= 0;	p1.ay = 0; 	p1.az = 0;
	p1.m  	= 0.56e-3;
	p1.vy 	= sqrt(G*(star.m+p1.m)/p1.x);
	p1.vx 	= 0;	p1.vz = 0;
	particles_add(p1); 
	
	struct particle p2;	// Planet 1
	p2.x 	= 1;	p2.y = 0; 	p2.z = 0;
	p2.ax 	= 0;	p2.ay = 0; 	p2.az = 0;
	p2.m  	= 1.89e-3;
	p2.vy 	= sqrt(G*(star.m+p2.m)/p2.x);
	p2.vx 	= 0;	p2.vz = 0;
	particles_add(p2); 

	tau_a = calloc(sizeof(double),N);
	tau_e = calloc(sizeof(double),N);

	tau_a[2] = 125663.;		// Migration timescale of planet 2 is 20000 years.
	tau_e[2] = tau_a[2]/100.; 	// Eccentricity damping timescale is 200 years (K=100). 
}

// Semi-major axis damping
void problem_adot(){
	for(int i=1;i<N;i++){
		if (tau_a[i]!=0){
			double tmpfac = dt/tau_a[i];
			// position
			struct particle* p = &(particles[i]);
			p->x  -= p->x*tmpfac;
			p->y  -= p->y*tmpfac;
			p->z  -= p->z*tmpfac;
			// velocity
			p->vx  += 0.5 * p->vx*tmpfac;
			p->vy  += 0.5 * p->vy*tmpfac;
			p->vz  += 0.5 * p->vz*tmpfac;
		}
	}
}

// Eccentricity damping
// This one is more complicated as it needs orbital elements to compute the forces.
void problem_edot(){
	for(int i=1;i<N;i++){
		if (tau_e[i]!=0){
			double d = dt/tau_e[i];
			struct particle* p = &(particles[i]);
			struct orbit o = tools_p2orbit(*p,particles[0].m);
			double rdot  = o.h/o.a/( 1. - o.e*o.e ) * o.e * sin(o.f);
			double rfdote = o.h/o.a/( 1. - o.e*o.e ) * ( 1. + o.e*cos(o.f) ) * (o.e + cos(o.f)) / (1.-o.e*o.e) / (1.+o.e*cos(o.f));
			//position
			double tmpfac = d * (  o.r/(o.a*(1.-o.e*o.e)) - (1.+o.e*o.e)/(1.-o.e*o.e));
			p->x -= tmpfac * p->x;
			p->y -= tmpfac * p->y;
			p->z -= tmpfac * p->z;
			//vx
			tmpfac = rdot/(o.e*(1.-o.e*o.e));
			p->vx -= d * o.e * (   cos(o.Omega) *      (tmpfac * cos(o.omega+o.f) - rfdote*sin(o.omega+o.f) )
						-cos(o.inc) * sin(o.Omega) * (tmpfac * sin(o.omega+o.f) + rfdote*cos(o.omega+o.f) ));
			//vy
			p->vy -= d * o.e * (   sin(o.Omega) *      (tmpfac * cos(o.omega+o.f) - rfdote*sin(o.omega+o.f) )
						+cos(o.inc) * cos(o.Omega) * (tmpfac * sin(o.omega+o.f) + rfdote*cos(o.omega+o.f) ));
			//vz
			p->vz -= d * o.e * (     sin(o.inc) *      (tmpfac * sin(o.omega+o.f) + rfdote*cos(o.omega+o.f) ));
		}
	}
}

void problem_inloop(){
}

void problem_output(){
	if (t>0){
		// The damping is done at the end of the timestep rather than in the loop
		// because we need the position and velocities at the same time to 
		// calculate orbital elements. We also update both the positions 
		// and velocities in these routines.
		problem_adot();
		problem_edot();
	}
	if(output_check(10000.*dt)){
		output_timing();
	}
	if(output_check(10.)){
		output_orbits_append("orbits.txt");
	}
}

void problem_finish(){
}
