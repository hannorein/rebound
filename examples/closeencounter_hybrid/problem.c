/**
 * Close Encounter with hybrid integrator (experimental)
 * 
 * This example integrates a densely packed planetary system 
 * which becomes unstable on a timescale of only a few orbits. 
 * This is a test case for the HYBRID integrator.
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

double e_init; // initial energy

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	r->dt = 0.012*2.*M_PI;						// initial timestep
	r->integrator = RB_IT_HYBRID;
	r->heartbeat  = heartbeat;

	struct reb_particle star;
	star.m = 1;
	star.x = 0; 	star.y = 0; 	star.z = 0;
	star.vx = 0; 	star.vy = 0; 	star.vz = 0;
	reb_add(r, star);
	
	// Add planets
	int N_planets = 3;
	for (int i=0;i<N_planets;i++){
		double a = 1.+.1*(double)i;		// semi major axis
		double v = sqrt(1./a); 					// velocity (circular orbit)
		struct reb_particle planet;
		planet.m = 2e-5; 
		planet.x = a; 	planet.y = 0; 	planet.z = 0;
		planet.vx = 0; 	planet.vy = v; 	planet.vz = 0;
		reb_add(r, planet); 
	}
	reb_move_to_com(r);				// This makes sure the planetary systems stays within the computational domain and doesn't drift.
	e_init = reb_tools_energy(r);
	system("rm -rf energy.txt");

	reb_integrate(r, INFINITY);
}

void heartbeat(struct reb_simulation* r){
	if (reb_output_check(r, 10.*2.*M_PI)){  
		reb_output_timing(r, 0);
	}
	if (reb_output_check(r, 2.*M_PI)){  
		FILE* f = fopen("energy.txt","a");
		reb_integrator_synchronize(r);
		double e = reb_tools_energy(r);
		fprintf(f,"%e %e\n",r->t, fabs((e-e_init)/e_init));
		fclose(f);
	}
}

