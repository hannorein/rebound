#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"

void problem_init(){
	// Setup particle structures
	init_particles(50);
	// Initial conditions
	for (int i =0;i<N;i++){
		particles[i].x = ((double)rand()/(double)RAND_MAX-0.5)*boxsize;
		particles[i].y = ((double)rand()/(double)RAND_MAX-0.5)*boxsize;
		particles[i].z = 0.1*((double)rand()/(double)RAND_MAX-0.5)*boxsize;
		particles[i].vx = 0;
		particles[i].vy = 0;
		particles[i].vz = 0;
		particles[i].ax = 0;
		particles[i].ay = 0;
		particles[i].az = 0;
		particles[i].m = 0.0001;
	}
	// Do not use any ghost boxes
	nghostx = 0;
	nghosty = 0;
	nghostz = 0;
}

void problem_output(){

}

void problem_finish(){

}
