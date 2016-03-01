/*
*  Stark Force in two body problem
*
*  Const. acceleration acting on planet at x direction
*  It would cause the eccentricity change between [0, 1]
*  by Meldonization. 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void stark_forces(struct reb_simulation* const r);
void heartbeat(struct reb_simulation* const r);
double tmax;

int main(int argc, char* argv[]) { 
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	r->dt 			= 1e-13;		// initial timestep.
	r->G			= 1;		// Gravitational constant
	r->integrator	= REB_INTEGRATOR_IAS15;

	// Setup callback function for velocity dependent forces.
	r->additional_forces 	= stark_forces;
	r->force_is_velocity_dependent = 0;
	// Setup callback function for outputs.
	r->heartbeat	= heartbeat;
	r->usleep		= 10000;		// Slow down integration (for visualization only)
	
	double mass_scale	= 1.;		// Some integrators have problems when changing the mass scale, IAS15 does not. 
	double size_scale	= 1.;		// Some integrators have problems when changing the size scale, IAS15 does not.
	
	struct reb_particle star = {0}; 
	star.m  = mass_scale;
	reb_add(r, star); 
	
	struct reb_particle p; 
	p.m  	= 0;	// massless
	p.x 	= size_scale;
	p.vy 	= sqrt(mass_scale/size_scale);
	reb_add(r, p); 

	reb_move_to_com(r);

	tmax	= 1e2*2.*M_PI;

	// Do the integration
	reb_integrate(r, tmax);

}


void stark_forces(struct reb_simulation* const r){
	double dragcoefficient = 0.01;
	struct reb_particle* const particles = r->particles;
	const int N = r->N;
	for (int i=1;i<N;i++){
		particles[i].ax += -dragcoefficient;
	}
}

void heartbeat(struct reb_simulation* const r){
		reb_output_timing(r, tmax);
}
