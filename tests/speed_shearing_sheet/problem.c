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


void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA = 1.;
	// Setup particle structures
	if (argc>1){
		boxsize = atof(argv[1]);
	}else{
		boxsize = 2;
	}
	root_nz 	= 4;
	dt 		= 1e-2;
	tmax 		= 10.*M_PI;
	int _N 		= round(50.*boxsize*boxsize);
	coefficient_of_restitution = 0.5;
	minimum_collision_velocity = 0.05;
	init_box();
	// Initial conditions
	for (int i =0;i<_N;i++){
		struct particle p;
		double vrand = 0.4*((double)rand()/(double)RAND_MAX-0.5);
		double phirand = 2.*M_PI*((double)rand()/(double)RAND_MAX-0.5);
		p.x = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_x;
		p.y = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
		p.z = .4*((double)rand()/(double)RAND_MAX-0.5);
		p.vx = 0;
		p.vy = -1.5*p.x*OMEGA+2.*vrand*cos(phirand);
		p.vz = vrand*sin(phirand);
		p.ax = 0;
		p.ay = 0;
		p.az = 0;
		p.m = 0.005;
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
	double error= 0;
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	double error_limit = 0.1;
	int pass = (error>error_limit?0:1);
	fprintf(of,"%d\t%e\t%e\t%e\t%d\t%e\t",N,dt,error,error_limit,pass,timing_final-timing_initial);
	fprintf(of,"N = %d",N);
	fprintf(of,"\n");
	fclose(of);
}
