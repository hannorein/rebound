#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	
	struct reb_particle p; 
	p.m  	= 1.;	
	reb_add(r, p); 

	double e = 1.e-14;

	struct reb_particle p1;
	p1.m = 0.;
	p1.x = 1.-e;
	p1.vy = 1.+e;
	reb_add(r,p1);

	struct reb_particle com = {0};
	struct reb_particle p2 = reb_tools_orbit_to_particle(r->G, com, 0., 2., 0., 0., 0., 0., 0.);
	reb_add(r,p2);
	//struct reb_orbit o1 = reb_tools_p2orbit(r->G, r->particles[1], r->particles[0]);
	//struct reb_orbit o2= reb_tools_p2orbit(r->G, r->particles[2], r->particles[0]);

	int i = reb_get_particle_index(&r->particles[1]);
	printf("index = %d\n", i);

}
