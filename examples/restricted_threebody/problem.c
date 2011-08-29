#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"

void problem_init(int argc, char* argv[]){
	// Setup constants
	boxsize 		= 8; 
	softening		= 1e-6;
	dt 			= 1.0e-2*2.*M_PI;
	N_active 		= 2; 	// Only the star and the planet have non-zero mass
	init_box();
	
	// Initial conditions for star
	struct particle star;
	star.x  = 0; 	star.y  = 0; 	star.z  = 0; 
	star.vx = 0; 	star.vy = 0; 	star.vz = 0;
	star.m  = 1;
	particles_add(star);

	// Initial conditions for planet
	double planet_e = 0.;
	struct particle planet;
	planet.x  = 1.-planet_e; 	planet.y  = 0; 				planet.z  = 0; 
	planet.vx = 0; 			planet.vy = sqrt(2./(1.-planet_e)-1.); 	planet.vz = 0;
	planet.m  = 1e-2;
	particles_add(planet);
	
	while(N<10000){
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
		particles_add(testparticle);
	}
}

void problem_output(){
	output_timing();
	if (output_check(2.*M_PI)){
		output_orbits("orbit.txt");
	}
}

void problem_finish(){
}

void problem_inloop(){
}

