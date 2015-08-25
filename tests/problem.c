/**
 * Velocity dependent drag force
 *
 * This is a very simple example on how to implement a velocity 
 * dependent drag force. The example uses the IAS15 integrator, which 
 * is ideally suited to handle non-conservative forces.
 * No gravitational forces or collisions are present.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "libreboundxf.h"
#include "tools.h"
#include "time.h"

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	
	struct reb_particle p; 
	p.m  	= 1.;	
	reb_add(r, p); 

	/*double e = 1.e-14;

	struct reb_particle p1;
	p1.m = 0.;
	p1.x = 1.-e;
	p1.vy = 1.+e;
	reb_add(r,p1);*/

	struct reb_particle p2 = reb_tools_init_orbit3d(r->G, 1., 0., 1., e, 0., 0., 0., 0.);
	reb_add(r,p2);
	//struct reb_orbit o1 = reb_tools_p2orbit(r->G, r->particles[1], r->particles[0]);
	struct reb_orbit o2= reb_tools_p2orbit(r->G, r->particles[2], r->particles[0]);
	
	/*printf("e = %.16e\n", o1.e);
	printf("w = %.16e\n", o1.omega);

	printf("e = %.16e\n", o2.e);
	printf("w = %.16e\n", o2.omega);

	printf("1-e*e = %.16f\n", 1.-e*e);

	double i = 0.1;
	double o = 0.3;
	double O = 0.2;
	struct reb_particle p3 = reb_tools_init_orbit3d(r->G, 1., 0., 1., 0., i, O, o, 0.);
	struct reb_particle p4 = reb_tools_init_orbit3d(r->G, 1., 0., 1., 0., -i, O, o, 0.);

	o1 = reb_tools_p2orbit(r->G, p3, r->particles[0]);
	o2 = reb_tools_p2orbit(r->G, p4, r->particles[0]);


	printf("%f\t%f\t%f\t%f\t%f\t%f\t\n", p3.x, p3.y, p3.z, p3.vx, p3.vy, p3.vz);
	printf("%f\t%f\t%f\t%f\t%f\t%f\t\n", p4.x, p4.y, p4.z, p4.vx, p4.vy, p4.vz);

	printf("%f\t%f\t%f\n", o1.inc, o1.Omega, o1.omega);
	printf("%f\t%f\t%f\n", o2.inc, o2.Omega, o2.omega);
	*/

}
