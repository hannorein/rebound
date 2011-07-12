#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "input.h"

// Example program using restart files.
// 
// First, run the program with './nbody'.
// Random initial conditions are created and
// a restart file is written once per orbit.
// Then, to restart the simulation, run the 
// program with './nbody --restart restart.bin'.
//

extern double OMEGA;
extern double coefficient_of_restitution;
extern double minimum_collision_velocity;

void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA = 1.;
	dt = 1e-4*M_PI*2.; 
	boxsize = 1;
	coefficient_of_restitution = 0.5;
	minimum_collision_velocity = .01;
	nghostx = 1;
	nghosty = 1;
	nghostz = 0;
	// Check if simulation is restarted
	if (input_check_restart(argc,argv)!=1){
		// Setup particle structures
		init_particles(50);
		// Initial conditions
		for (int i =0;i<N;i++){
			particles[i].x  = ((double)rand()/(double)RAND_MAX-0.5)*boxsize;
			particles[i].y  = ((double)rand()/(double)RAND_MAX-0.5)*boxsize;
			particles[i].z  = 0.1*((double)rand()/(double)RAND_MAX-0.5)*boxsize;
			particles[i].vx = 0;
			particles[i].vy = -1.5*particles[i].x*OMEGA;
			particles[i].vz = 0;
			particles[i].ax = 0;
			particles[i].ay = 0;
			particles[i].az = 0;
			particles[i].m  = 0.01;
			particles[i].r  = 0.05;
		}
	}
}

void problem_inloop(){

}

void problem_output(){
	// Override the restartfile every orbit.
	if (output_check(1e-0*2.*M_PI&&t!=0)){
		output_binary("restart.bin");
	}
}

void problem_finish(){

}
