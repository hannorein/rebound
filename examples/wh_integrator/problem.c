#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"

// The Wisdom Holman integrator assumes a fixed mass
// at r=0 with mass m=1.
// 

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt = 1e-1;
	boxsize = 3;
	// Setup particle structures
	init_particles(1);
	// Initial conditions
	particles[0].x  = 1;
	particles[0].y  = 0;
	particles[0].z  = 0;
	particles[0].vx = 0;
	particles[0].vy = 1;
	particles[0].vz = 0;
	particles[0].ax = 0;
	particles[0].ay = 0;
	particles[0].az = 0;
	particles[0].m  = 0;
	// Do not use any ghost boxes
	nghostx = 0;
	nghosty = 0;
	nghostz = 0;
}

void problem_inloop(){

}

void problem_output(){

}

void problem_finish(){

}
