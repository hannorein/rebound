#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "collisions.h"

extern double OMEGA;
extern double OMEGAZ;
extern double coefficient_of_restitution; 
extern double minimum_collision_velocity;

void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA 				= 1.;
	OMEGAZ 				= 3.6;
	dt				= 1e-2*2.*M_PI/OMEGA;
	double particle_r 		= 1;
	double tau			= 1;
	coefficient_of_restitution 	= 0.5;
//	minimum_collision_velocity 	= 0.001*particle_r*OMEGA;
	boxsize 			= 10;
	root_nx = 100; 	root_ny = 1; 	root_nz = 2;
	nghostx = 1; 	nghosty = 1; 	nghostz = 0;
	init_box();
	// Initial conditions
	double _N = tau * boxsize_x * boxsize_y/(M_PI*particle_r *particle_r);
	while (N<_N){
		struct particle p;
		p.x 	= ((double)rand()/(double)RAND_MAX-0.5)*boxsize_x;
		p.y 	= ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
		p.z 	= 4.0*((double)rand()/(double)RAND_MAX-0.5)*particle_r;
		p.vx 	= 0;
		p.vy 	= -1.5*p.x*OMEGA;
		p.vz 	= 0;
		p.ax 	= 0;
		p.ay 	= 0;
		p.az 	= 0;
		p.m 	= 0.0001;
		p.r 	= particle_r;
		particles_add(p);
	}
}

void problem_inloop(){

}

void problem_output(){
	if (output_check(dt)){
		output_timing();
	}
}

void problem_finish(){

}
