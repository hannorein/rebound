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

extern double opening_angle2;

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
	double surfacedensity 	= 0.06; 
	double particle_mass 	= 5e-6;
	int _N = round(surfacedensity*M_PI*(boxsize*boxsize/4./sqrt(1.2)-boxsize*boxsize/400.)/particle_mass);
	// Initial conditions
	struct particle star;
	star.x 		= 0; star.y 	= 0; star.z	= 0;
	star.vx 	= 0; star.vy 	= 0; star.vz 	= 0;
	star.ax 	= 0; star.ay 	= 0; star.az 	= 0;
	star.m 		= 1;
	particles_add(star);
	while(N<_N){
		struct particle pt;
		pt.x 		= ((double)rand()/(double)RAND_MAX-0.5)*boxsize_x;
		pt.y 		= ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
		pt.z 		= ((double)rand()/(double)RAND_MAX-0.5)*boxsize_z*0.1;
		double a	 = sqrt(pt.x*pt.x + pt.y*pt.y);
		if (a>boxsize/2./1.2 || a<boxsize/20.) continue;
		double phi 	= atan2(pt.y,pt.x);
		double vkep 	= sqrt(G*star.m/(a));
		pt.vx 		=  vkep * sin(phi);
		pt.vy 		= -vkep * cos(phi);
		pt.vz 		= 0;
		pt.ax 		= 0;
		pt.ay 		= 0;
		pt.az 		= 0;
		pt.m 		= particle_mass;
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
