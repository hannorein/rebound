/**
 * A very simple test problem
 * 
 * We first create a REBOUND simulation, then we add 
 * two particles and integrate the system for 100 time 
 * units.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

void heartbeat(struct reb_simulation* r){
	printf("%f\n",r->t);
}

int main(int argc, char* argv[]) {
	struct reb_simulation* r = reb_create_simulation();
	r->dt = 0.1;
	r->heartbeat = heartbeat;

	struct reb_particle p1;
	p1.x = 0;  p1.y = 0;  p1.z = 0; 
	p1.vx = 0; p1.vy = 0; p1.vz = 0; 
	p1.m = 1.;
	reb_add(r, p1);
	
	struct reb_particle p2;
	p2.x = 1;  p2.y = 0;  p2.z = 0; 
	p2.vx = 0; p2.vy = 1; p2.vz = 0; 
	p2.m = 0.;
	reb_add(r, p2);

	reb_integrate(r,100.);
}

