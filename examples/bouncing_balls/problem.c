/**
 * Bouncing balls.
 * 
 * This example is a simple test of collision detection 
 * methods.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	r->integrator	= REB_INTEGRATOR_LEAPFROG;
	r->gravity	= REB_GRAVITY_BASIC;
	r->collision	= REB_COLLISION_DIRECT;
	r->dt = 1e-2;

	reb_configure_box(r, 3.0, 1, 1, 1);
	
	// Initial conditions
	{
		struct reb_particle p;
		p.x  = 1; p.y  = 1; p.z  = 1;
		p.vx = 0; p.vy = 0; p.vz = 0;
		p.ax = 0; p.ay = 0; p.az = 0;
		p.m  = 1;
		p.r  = 0.1;
		reb_add(r, p);
	}
	{
		struct reb_particle p;
		p.x  = -1; p.y  = -1; p.z  = -1;
		p.vx =  0; p.vy =  0; p.vz =  0;
		p.ax =  0; p.ay =  0; p.az =  0;
		p.m  = 1;
		p.r  = 0.1;
		reb_add(r, p);
	}

	reb_integrate(r, INFINITY);
}

