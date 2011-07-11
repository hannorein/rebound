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
	boxsize = 3;
	coefficient_of_restitution = 1; // elastic collisions

	// Setup particle structures
	init_particles(4);
	// Initial conditions
	particles[0].x  = 1;
	particles[0].y  = 1;
	particles[0].z  = 0;
	particles[0].vx = 0;
	particles[0].vy = 0;
	particles[0].vz = 0;
	particles[0].ax = 0;
	particles[0].ay = 0;
	particles[0].az = 0;
	particles[0].m  = 1;
	particles[0].r  = 0.1;
	particles[1].x  = -1;
	particles[1].y  = -1;
	particles[1].z  = 0;
	particles[1].vx = 0;
	particles[1].vy = 0;
	particles[1].vz = 0;
	particles[1].ax = 0;
	particles[1].ay = 0;
	particles[1].az = 0;
	particles[1].m  = 1;
	particles[1].r  = 0.1;
	particles[2].x  = 1;
	particles[2].y  = -1;
	particles[2].z  = 0;
	particles[2].vx = 0;
	particles[2].vy = 0;
	particles[2].vz = 0;
	particles[2].ax = 0;
	particles[2].ay = 0;
	particles[2].az = 0;
	particles[2].m  = 1;
	particles[2].r  = 0.1;
	particles[3].x  = -1;
	particles[3].y  = 1;
	particles[3].z  = 0;
	particles[3].vx = 0;
	particles[3].vy = 0;
	particles[3].vz = 0;
	particles[3].ax = 0;
	particles[3].ay = 0;
	particles[3].az = 0;
	particles[3].m  = 1;
	particles[3].r  = 0.1;
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
