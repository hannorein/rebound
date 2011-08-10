#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"

double e =0;
void problem_init(int argc, char* argv[]){
	if (argc>1){
		e = atof(argv[1]);
	}
	// Setup constants
	dt 		= 1e-3*2.*M_PI;
	tmax 		= 2.*M_PI;
	boxsize 	= 3e9;
	N_active_first 	= 1;
	N_active_last 	= 2;
	init_box();
	// Initial conditions
	struct particle p0;
	p0.x  = 0;
	p0.y  = 0;
	p0.z  = 0;
	p0.vx = 0;
	p0.vy = 0;
	p0.vz = 0;
	p0.ax = 0;
	p0.ay = 0;
	p0.az = 0;
	p0.m  = 1;
	particles_add(p0);
	struct particle p1;
	p1.x  = (1.-e);
	p1.y  = 0;
	p1.z  = 0;
	p1.vx = 0;
	p1.vy = sqrt((1.+e)/(1.-e));
	p1.vz = 0;
	p1.ax = 0;
	p1.ay = 0;
	p1.az = 0;
	p1.m  = 0;
	particles_add(p1);
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
	FILE* of = fopen("error.txt","a+"); 
	struct particle p = particles[1];
	double dx = fabs(p.x-(1.-e));
	double dy = fabs(p.y-0.);
	double dz = fabs(p.z-0.);
	double error= sqrt(dx*dx+dy*dy+dz*dz);
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	double error_limit = 5e-11;
	int pass = (error>error_limit?0:1);
	fprintf(of,"%d\t%e\t%e\t%e\t%d\t%e\t",N,dt,error,error_limit,pass,timing_final-timing_initial);
	fprintf(of,"eccentricity = %e",e);
	fprintf(of,"\n");
	fclose(of);
}
