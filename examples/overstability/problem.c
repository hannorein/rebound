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
	OMEGA = 1.;
	OMEGAZ = 3.6;
	// Setup particle structures
	boxsize = 30;
	init_particles(1000);
	dt = 2e-2*2.*M_PI;
	double particle_size = 0.2;
	coefficient_of_restitution = 0.5;
	minimum_collision_velocity = 0.01*particle_size*OMEGA;
	collisions_max_r = particle_size;
	// Initial conditions
	for (int i =0;i<N;i++){
		particles[i].x = ((double)rand()/(double)RAND_MAX-0.5)*boxsize;
		particles[i].y = ((double)rand()/(double)RAND_MAX-0.5)*boxsize;
		particles[i].z = 0.01*((double)rand()/(double)RAND_MAX-0.5)*boxsize;
		particles[i].vx = 0;
		particles[i].vy = -1.5*particles[i].x*OMEGA;
		particles[i].vz = 0;
		particles[i].ax = 0;
		particles[i].ay = 0;
		particles[i].az = 0;
		particles[i].m = 0.0001;
		particles[i].r = particle_size;
	}
	// Do use ghost boxes in x and y
	nghostx = 1;
	nghosty = 1;
	nghostz = 0;
}

void problem_inloop(){

}

void problem_output(){
	if (output_check(1e1*2.*M_PI)){
		output_timing();
	}
}

void problem_finish(){

}
