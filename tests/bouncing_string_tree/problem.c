#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"

extern double coefficient_of_restitution; 

int loops = 1;
int loopsleft;

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt = 1e-4;
	tmax = 1;
	if (argc>1){
		loops = atoi(argv[1]);
	}
	loops*=2;
	loopsleft = loops;

	tmax = 1.2*(double)loops;
	boxsize = 3;
	root_nx = 10;
	root_ny = 3;
	root_nz = 1;
	coefficient_of_restitution = 1; // elastic collisions

	// Setup particle structures
	init_box();
	// Initial conditions
	for (int i=0;i<10;i++){
		struct particle pt;
		pt.x  = -boxsize_x/2.+boxsize_x*(double)i/10.;
		pt.y  = 0;
		pt.z  = 0;
		pt.vx = (i==0?10:0);
		pt.vy = 0;
		pt.vz = 0;
		pt.ax = 0;
		pt.ay = 0;
		pt.az = 0;
		pt.m  = 1;
		pt.r  = 1;
		particles_add(pt);
	}
	// Do not use any ghost boxes
	nghostx = 1;
	nghosty = 1;
	nghostz = 0;
	printf("init done\n");
}

double tcontact =0 ;
double vlast = 0;
extern double collisions_plog;
void problem_inloop(){
	if (particles[0].vx!=vlast){
		vlast = particles[0].vx;
		if (vlast>0){
			tcontact = t;
		}
	}
}

void problem_output(){

}

void problem_finish(){
	FILE* of = fopen("error.txt","a+"); 
	double error= fabs(1.1*(double)loops-tcontact);
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	double error_limit = 2e-3;
	int pass = (error>error_limit?0:1);
	fprintf(of,"%d\t%e\t%e\t%e\t%d\t%e\t",N,dt,error,error_limit,pass,timing_final-timing_initial);
	fprintf(of,"loops = %d",loops/2);
	fprintf(of,"\n");
	fclose(of);
}
