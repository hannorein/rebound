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
#include "tools.h"


double total_energy();
double energy_initial;

void problem_init(int argc, char* argv[]){
	boxsize 	= 20;
	dt 		= 1e-3;
	tmax		= 10;
	init_box();
	// Initial conditions
	while(N<100){
		struct particle p;
		p.x  = tools_uniform(-boxsize_x/2.,boxsize/2.);
		p.y  = tools_uniform(-boxsize_x/2.,boxsize/2.);
		p.z  = tools_uniform(-boxsize_x/2.,boxsize/2.);
		p.vx = tools_normal(1);
		p.vy = tools_normal(1);
		p.vz = tools_normal(1);
		p.ax = 0;
		p.ay = 0;
		p.az = 0;
		p.m  = 1;
		p.r  = 1;
		particles_add(p);
	}
	energy_initial = total_energy();
	// No ghost boxes 
	nghostx = 1;
	nghosty = 1;
	nghostz = 1;
}

double total_energy(){
	double energy = 0;
	for (int i=0;i<N;i++){
		struct particle p1 = particles[i];
		energy += 0.5*p1.m*(p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz);
	}
	return energy;
}

void problem_inloop(){
}

void problem_output(){
	output_timing();
}

void problem_finish(){
	FILE* of = fopen("error.txt","a+"); 
	double error= fabs((energy_initial-total_energy())/energy_initial);
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	double error_limit = 1e-10;
	int pass = (error>error_limit?0:1);
	fprintf(of,"%d\t%e\t%e\t%e\t%d\t%e\t",N,dt,error,error_limit,pass,timing_final-timing_initial);
	fprintf(of,"N = %d",N);
	fprintf(of,"\n");
	fclose(of);
}
