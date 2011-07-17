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
	dt = 1e-3*2.*M_PI;
	tmax = 2.*M_PI;
	boxsize = 3e9;
	// Setup particle structures
	init_particles(2);
	N_active_first = 1;
	N_active_last = 2;
	// Initial conditions
	particles[0].x  = 0;
	particles[0].y  = 0;
	particles[0].z  = 0;
	particles[0].vx = 0;
	particles[0].vy = 0;
	particles[0].vz = 0;
	particles[0].ax = 0;
	particles[0].ay = 0;
	particles[0].az = 0;
	particles[0].m  = 1;
	particles[1].x  = (1.-e);
	particles[1].y  = 0;
	particles[1].z  = 0;
	particles[1].vx = 0;
	particles[1].vy = sqrt((1.+e)/(1.-e));
	particles[1].vz = 0;
	particles[1].ax = 0;
	particles[1].ay = 0;
	particles[1].az = 0;
	particles[1].m  = 0;
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
