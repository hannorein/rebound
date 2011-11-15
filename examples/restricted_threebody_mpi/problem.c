#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "communication_mpi.h"

void problem_init(int argc, char* argv[]){
	// Setup constants
	boxsize 		= 8; 
	softening		= 1e-6;
	dt 			= 1.0e-2*2.*M_PI;
	N_active 		= 2; 	// Only the star and the planet have non-zero mass
	root_nx	= 2; root_ny	= 2; root_nz	= 1; 
	init_box();
	
	// Initial conditions for star
	struct particle star;
	star.x  = 0; 	star.y  = 0; 	star.z  = 0; 
	star.vx = 0; 	star.vy = 0; 	star.vz = 0;
	star.m  = 1;

	// Initial conditions for planet
	double planet_e = 0.;
	struct particle planet;
	planet.x  = 1.-planet_e; 	planet.y  = 0; 				planet.z  = 0; 
	planet.vx = 0; 			planet.vy = sqrt(2./(1.-planet_e)-1.); 	planet.vz = 0;
	planet.m  = 1e-2;
	
	int _N = 10000;

	// Move to centre of mass frame
	double com_x  = (star.x*star.m  + planet.x*planet.m) /(star.m+planet.m);
	double com_y  = (star.y*star.m  + planet.y*planet.m) /(star.m+planet.m);
	double com_z  = (star.z*star.m  + planet.z*planet.m) /(star.m+planet.m);
	double com_vx = (star.vx*star.m + planet.vx*planet.m)/(star.m+planet.m);
	double com_vy = (star.vy*star.m + planet.vy*planet.m)/(star.m+planet.m);
	double com_vz = (star.vz*star.m + planet.vz*planet.m)/(star.m+planet.m);
	planet.x  -= com_x; 	planet.y  -= com_y; 	planet.z  -= com_z;
	planet.vx -= com_vx; 	planet.vy -= com_vy; 	planet.vz -= com_vz;
	star.x    -= com_x; 	star.y    -= com_y; 	star.z    -= com_z;
	star.vx   -= com_vx; 	star.vy   -= com_vy; 	star.vz   -= com_vz;
	// Add active particles on all nodes
	particles_add(star);
	particles_add(planet);
#ifdef MPI
	_N /= mpi_num;
#endif	// MPI
	
	while(N<_N){
		double x 	= ((double)rand()/(double)RAND_MAX-0.5)*boxsize*0.9;
		double y 	= ((double)rand()/(double)RAND_MAX-0.5)*boxsize*0.9;
		double a 	= sqrt(x*x+y*y);
		double phi 	= atan2(y,x);
		if (a<.1) continue;
		if (a>boxsize_x/2.*0.9) continue;

		double vkep = sqrt(G*star.m/a);
		struct particle testparticle;
		testparticle.x  = x;
		testparticle.y  = y; 
		testparticle.z  = 1.0e-2*x*((double)rand()/(double)RAND_MAX-0.5);
		testparticle.vx = -vkep*sin(phi);
		testparticle.vy = vkep*cos(phi);
		testparticle.vz = 0;
		testparticle.ax = 0; 
		testparticle.ay = 0; 
		testparticle.az = 0;
		testparticle.m  = 0;

		// Only add local particles. This is quick and dirty. 
		// Definitely not the most efficient way.
#ifdef MPI
		int rootbox = particles_get_rootbox_for_particle(testparticle);
		int root_n_per_node = root_n/mpi_num;
		int proc_id = rootbox/root_n_per_node;
		if (proc_id != mpi_id) continue;
#endif	//MPI
		
		particles_add(testparticle);
	}
}

void problem_output(){
	if (output_check(2.*M_PI)){
		output_timing();
	}
	if (output_check(2.*M_PI)){
		output_ascii("positions.txt");
	}
}

void problem_finish(){
}

void problem_inloop(){
}

