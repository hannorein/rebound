#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "communication_mpi.h"
#include "tree.h"
#include "tools.h"

extern double opening_angle2;
extern int WH_SELFGRAVITY_ENABLED;
extern int Nmax;

void problem_init(int argc, char* argv[]){
	// Setup constants
	opening_angle2	= 1.5;
	G 		= 1;		
	softening 	= 0.01;		
	dt 		= 3e-3;
	boxsize 	= 1.2;
	root_nx = 1; root_ny = 1; root_nz = 1;
	nghostx = 0; nghosty = 0; nghostz = 0; 		
	init_box();
	// Setup particles
	double disc_mass = 2e-1;
	int _N = 10000;
	// Initial conditions
	struct particle star;
	star.x 		= 0; star.y 	= 0; star.z	= 0;
	star.vx 	= 0; star.vy 	= 0; star.vz 	= 0;
	star.ax 	= 0; star.ay 	= 0; star.az 	= 0;
	star.m 		= 1;
#ifdef INTEGRATOR_WH
	// Insert particle manually. Don't add it to tree.
	WH_SELFGRAVITY_ENABLED 	= 1;
	N_active_first 		= 1;
	Nmax 			+= 128;
	particles 		= realloc(particles,sizeof(struct particle)*Nmax);
	particles[N] 		= star;
	N++;
#else // INTEGRATOR_WH
	particles_add(star);
#endif // INTEGRATOR_WH
	while(N<_N){
		struct particle pt;
		double a	= tools_powerlaw(boxsize/10.,boxsize/2./1.2,-1.5);
		double phi 	= tools_uniform(0,2.*M_PI);
		pt.x 		= a*cos(phi);
		pt.y 		= a*sin(phi);
		pt.z 		= a*tools_normal(0.001);
		double mu 	= star.m + disc_mass * (pow(a,-3./2.)-pow(boxsize/10.,-3./2.))/(pow(boxsize/2./1.2,-3./2.)-pow(boxsize/10.,-3./2.));
		printf("a = %f\tmu = %f\n",a,mu);
		double vkep 	= sqrt(G*mu/a);
		pt.vx 		=  vkep * sin(phi);
		pt.vy 		= -vkep * cos(phi);
		pt.vz 		= 0;
		pt.ax 		= 0;
		pt.ay 		= 0;
		pt.az 		= 0;
		pt.m 		= disc_mass/(double)_N;
		particles_add(pt);
	}
}

void problem_inloop(){
}

void problem_output(){
	if (output_check(10.0*dt)) output_timing();
}

void problem_finish(){
}
