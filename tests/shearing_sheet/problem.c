#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"

extern double OMEGA;
extern double coefficient_of_restitution;
extern double minimum_collision_velocity;
int N_init = 10;

void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA 		= 1.;
	// Setup particle structures
	boxsize 	= 1;
	root_nx		= 1;
	root_ny 	= 4;
	root_nz 	= 1;
	coefficient_of_restitution = 0.5;
	minimum_collision_velocity = 0.005;
	dt 		= 1e-2;
	tmax		= 10.*2.*M_PI;
	init_box();
	// Initial conditions
	for (int i =0;i<N_init;i++){
		struct particle p;
		p.x = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_x;
		p.y = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
		p.z = 0;//0.1*((double)rand()/(double)RAND_MAX-0.5)*boxsize;
		p.vx = 0;
		p.vy = -1.5*p.x*OMEGA;
		p.vz = 0;
		p.ax = 0;
		p.ay = 0;
		p.az = 0;
		p.m = 0.0001;
		p.r = 0.1;
		particles_add(p);
	}
	// Do use ghost boxes in x and y
	nghostx = 1;
	nghosty = 1;
	nghostz = 0;
}

void problem_inloop(){

}

void problem_output(){
}

void problem_finish(){
	FILE* of = fopen("error.txt","a+"); 
	double error= N_init-N;
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	double error_limit = 0.1;
	int pass = (error>error_limit?0:1);
	fprintf(of,"%d\t%e\t%e\t%e\t%d\t%e\t",N,dt,error,error_limit,pass,timing_final-timing_initial);
	fprintf(of,"N_init = %d",N_init);
	fprintf(of,"\n");
	fclose(of);

}
