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
#include "tools.h"
#include "output.h"
#include "particle.h"
#include "boundaries.h"

double* tau_a; 	/**< Migration timescale in years for all particles */
double* tau_e; 	/**< Eccentricity damping timescale in years for all particles */

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 		= 1e-3;
	boxsize 	= 3;
	init_box();

	// Initial conditions
	struct particle star;
	star.x  = 0; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = 0; star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = 0.32;
	particles_add(star); 
	
	struct particle p1;
	p1.x 	= 0.5;	p1.y = 0;	p1.z = 0;
	p1.ax 	= 0;	p1.ay = 0; 	p1.az = 0;
	p1.m  	= 0.56e-3;
	p1.vy 	= sqrt(G*(star.m+p1.m)/p1.x);
	p1.vx 	= 0;	p1.vz = 0;
	particles_add(p1); 
	
	struct particle p2;
	p2.x 	= 1;	p2.y = 0; 	p2.z = 0;
	p2.ax 	= 0;	p2.ay = 0; 	p2.az = 0;
	p2.m  	= 1.89e-3;
	p2.vy 	= sqrt(G*(star.m+p2.m)/p2.x);
	p2.vx 	= 0;	p2.vz = 0;
	particles_add(p2); 

	tau_a = calloc(sizeof(double),N);
	tau_e = calloc(sizeof(double),N);

	//tau_a[1] = 1000;
	//tau_a[2] = 1000;
	//tau_e[2] = 1000;
}

void problem_adot(){
	for(int i=0;i<N;i++){
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


void problem_edot(){
	for(int i=0;i<N;i++){
		if (tau_e[i]!=0){
			struct particle* p = &(particles[i]);
			struct orbit o = tools_p2orbit(*p,particles[0].m);
			double rdot  = o.h/o.a/( 1. - o.e*o.e ) * o.e * sin(o.f);
			double rfdote = o.h/o.a/( 1. - o.e*o.e ) * ( 1. + o.e*cos(o.f) ) * (o.e + cos(o.f)) / (1.-o.e*o.e) / (1.+o.e*cos(o.f));
			//position
			double tmpfac =   o.r/(o.a*(1.-o.e*o.e)) - (1.+o.e*o.e)/(1.-o.e*o.e);
			p->x -= tmpfac * p->x;
			p->y -= tmpfac * p->y;
			p->z -= tmpfac * p->z;
			//vx
			tmpfac = rdot/(o.e*(1.-o.e*o.e));
			p->vx -= o.e * ( cos(o.Omega) *      (tmpfac * cos(o.omega+o.f) - rfdote*sin(o.omega+o.f) )
					-cos(o.inc) * sin(o.Omega) * (tmpfac * sin(o.omega+o.f) + rfdote*cos(o.omega+o.f) ));
			//vy
			p->vy -= o.e * ( sin(o.Omega) *      (tmpfac * cos(o.omega+o.f) - rfdote*sin(o.omega+o.f) )
					+cos(o.inc) * cos(o.Omega) * (tmpfac * sin(o.omega+o.f) + rfdote*cos(o.omega+o.f) ));
			//vz
			p->vz -= o.e * ( sin(o.inc) *        (tmpfac * sin(o.omega+o.f) + rfdote*cos(o.omega+o.f) ));
		}
	}
}

void problem_inloop(){
	if(output_check(100.*dt)){
		output_timing();
	}
}

void problem_output(){
	if (t>0){
		problem_adot();
		problem_edot();
	}
}

void problem_finish(){
}
