/**
 * @file 	problem.c
 * @brief 	Example problem: Velocity dependent drag force
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This is a very simple example on how to implement a velocity 
 * dependent drag force. The example uses the IAS15 integrator, which 
 * is ideally suited to handle non-conservative forces.
 * No gravitational forces or collisions are present.
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
#include "rebound.h"
#include "tools.h"
#include "output.h"
#include "particle.h"
#include "integrator.h"

void additional_forces(struct reb_simulation* const r);
void heartbeat(struct reb_simulation* const r);

double tmax = 40.;

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	r->dt 			= 1e-4;		// initial timestep.
	r->integrator		= RB_IT_IAS15;
	r->gravity		= RB_GT_NONE;

	// Setup callback function for velocity dependent forces.
	r->additional_forces 	= additional_forces;
	r->force_is_velocitydependent = 1;
	// Setup callback function for outputs.
	r->heartbeat		= heartbeat;
	
	struct reb_particle p; 
	p.m  = 0;	// massless
	p.x = 1; 	p.y = 0; 	p.z = 0;
	p.vx = -1; 	p.vy = 0; 	p.vz = 0;
	p.ax = 0; 	p.ay = 0; 	p.az = 0;
	reb_add(r, p); 

	// Delete previous output
	system("rm -v r.txt");	

	// Do the integration
	reb_integrate(r, tmax);
}

void additional_forces(struct reb_simulation* const r){
	// Simplest velocity dependent drag force.
	double dragcoefficient = 1;
	struct reb_particle* const particles = r->particles;
	const int N = r->N;
	for (int i=0;i<N;i++){
		particles[i].ax = -dragcoefficient*particles[i].vx;
		particles[i].ay = -dragcoefficient*particles[i].vy;
		particles[i].az = -dragcoefficient*particles[i].vz;
	}
}

void heartbeat(struct reb_simulation* const r){
	// Output some information to the screen every 100th timestep
	if(reb_output_check(r, 100.*r->dt)){
		reb_output_timing(r, tmax);
	}
	// Output the particle position to a file every timestep.
	const struct reb_particle* const particles = r->particles;
	FILE* f = fopen("r.txt","a");
	fprintf(f,"%e\t%e\t%e\n",r->t,particles[0].x, particles[1].vx);
	fclose(f);
}
