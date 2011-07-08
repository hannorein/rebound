#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"

extern double OMEGA;

void problem_init(){
	// Setup constants
	OMEGA = 1.;
	// Setup particle structures
	init_particles(50);
	// Initial conditions
	for (int i =0;i<N;i++){
		particles[i].x = ((double)rand()/(double)RAND_MAX-0.5)*boxsize;
		particles[i].y = ((double)rand()/(double)RAND_MAX-0.5)*boxsize;
		particles[i].z = 0.1*((double)rand()/(double)RAND_MAX-0.5)*boxsize;
		particles[i].vx = 0;
		particles[i].vy = -1.5*particles[i].x*OMEGA;
		particles[i].vz = 0;
		particles[i].ax = 0;
		particles[i].ay = 0;
		particles[i].az = 0;
		particles[i].m = 0.0001;
	}
	// Do use ghost boxes in x and y
	nghostx = 1;
	nghosty = 1;
	nghostz = 0;
}

void problem_output(){

}

void problem_finish(){

}
