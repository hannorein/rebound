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
	dt 	= 1e-3;
	boxsize = 3;
	coefficient_of_restitution = 1; // elastic collisions
	
	int problem_id = 1;
	if (argc>1){			
		problem_id = atoi(argv[1]);
	}

	// Setup particle structures
	init_box();
	// Initial conditions
	struct particle p;

	switch (problem_id){
		case 1:
			nghostx = 1; nghosty = 1; nghostz = 0;
			p.x  = 1; p.y  = 1; p.z  = 0;
			p.vx = 0; p.vy = 0; p.vz = 0;
			p.ax = 0; p.ay = 0; p.az = 0;
			p.m  = 1;
			p.r  = 0.1;
			particles_add(p);
			p.x  = -1; p.y  = -1; p.z  = 0;
			p.vx = 0;  p.vy = 0;  p.vz = 0;
			p.ax = 0;  p.ay = 0;  p.az = 0;
			p.m  = 1;
			p.r  = 0.1;
			particles_add(p);
			p.x  = 1; p.y  = -1; p.z  = 0;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 1;
			p.r  = 0.1;
			particles_add(p);
			p.x  = -1; p.y  = 1; p.z  = 0;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 1;
			p.r  = 0.1;
			particles_add(p);
			break;
		case 2:
			nghostx = 0; nghosty = 0; nghostz = 0;
			p.x  = 0; p.y  = 0; p.z  = 0;
			p.vx = 0; p.vy = 0; p.vz = 0;
			p.ax = 0; p.ay = 0; p.az = 0;
			p.m  = 1;
			p.r  = 0.5;
			particles_add(p);
			p.x  = 1; p.y  = 1; p.z  = 0;
			p.vx = 0; p.vy = 0; p.vz = 0;
			p.ax = 0; p.ay = 0; p.az = 0;
			p.m  = 0.008;
			p.r  = 0.1;
			particles_add(p);
			p.x  = -1; p.y  = -1; p.z  = 0;
			p.vx = 0;  p.vy = 0;  p.vz = 0;
			p.ax = 0;  p.ay = 0;  p.az = 0;
			p.m  = 0.008;
			p.r  = 0.1;
			particles_add(p);
			p.x  = 1; p.y  = -1; p.z  = 0;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 0.008;
			p.r  = 0.1;
			particles_add(p);
			p.x  = -1; p.y  = 1; p.z  = 0;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 0.008;
			p.r  = 0.1;
			particles_add(p);
			p.x  = 0;  p.y  = 0; p.z  = 1.3;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 0.008;
			p.r  = 0.1;
			particles_add(p);
			p.x  = 0;  p.y  = 0; p.z  =-1.3;
			p.vx = 0;  p.vy = 0; p.vz = 0;
			p.ax = 0;  p.ay = 0; p.az = 0;
			p.m  = 0.008;
			p.r  = 0.1;
			particles_add(p);
			break;


	}
	printf("init done\n");
}

void problem_inloop(){

}

void problem_output(){

}

void problem_finish(){

}
