#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"

extern double coefficient_of_restitution; 
void problem_init(int argc, char* argv[]){
	// Setup constants
	dt = 1e-3;
	tmax = 10000;
	boxsize_x = 30;
	boxsize_y = 10;
	boxsize_z = 3;
	coefficient_of_restitution = 1; // elastic collisions

	// Setup particle structures
	init_particles(10);
	// Initial conditions
	for (int i=0;i<N;i++){
		particles[i].x  = -boxsize_x/2.+boxsize_x*(double)i/(double)N;
		particles[i].y  = 0;
		particles[i].z  = 0;
		particles[i].vx = 0;
		particles[i].vy = 0;
		particles[i].vz = 0;
		particles[i].ax = 0;
		particles[i].ay = 0;
		particles[i].az = 0;
		particles[i].m  = 1;
		particles[i].r  = 1;
	}
	particles[0].vx = 20;
	// Do not use any ghost boxes
	nghostx = 1;
	nghosty = 1;
	nghostz = 0;
	printf("init done\n");
}

void problem_inloop(){

}

void problem_output(){

}

void problem_finish(){

}
