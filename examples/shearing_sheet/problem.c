#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"

extern double OMEGA;
extern double coefficient_of_restitution;
extern double minimum_collision_velocity;

void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA = 1.;
	// Setup particle structures
	boxsize_x = 2;
	boxsize_y = 4;
	boxsize_z = 1;
	init_particles(50);
	coefficient_of_restitution = 0.5;
	minimum_collision_velocity = 0.005;
	dt = 1e-3;
	// Initial conditions
	for (int i =0;i<N;i++){
		particles[i].x = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_x;
		particles[i].y = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
		particles[i].z = 0;//0.1*((double)rand()/(double)RAND_MAX-0.5)*boxsize;
		particles[i].vx = 0;
		particles[i].vy = -1.5*particles[i].x*OMEGA;
		particles[i].vz = 0;
		particles[i].ax = 0;
		particles[i].ay = 0;
		particles[i].az = 0;
		particles[i].m = 0.0001;
		particles[i].r = 0.1;
	}
	// Do use ghost boxes in x and y
	nghostx = 1;
	nghosty = 1;
	nghostz = 0;
}

void problem_inloop(){

}

void problem_output(){
	if (output_check(1e-1*2.*M_PI)){
		output_append_velocity_dispersion("veldisp.txt");
	}
}

void problem_finish(){

}
