#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"

extern double OMEGA;
double cutoff_radius;
double g = 1;
double amplitude_start = 1e-10;
double amplitude_finish = 1e-5;

void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA = 1.;
	boxsize = 4; // = Number of particles
	dt = 1.12342535151342e-3*2.*M_PI;
	softening = 0;
	tmax = 100.*2.*M_PI;
	if (argc>=2){
		g = atof(argv[1]);
	}else{
		fprintf(stderr,"No g given!\n");
		exit(-1);
	}
	// Setup particle structures
	init_particles((int)boxsize);
	// Initial conditions
	int i=0;
	for (double y =-boxsize/2.0+boxsize/(double)N*0.5;i<N;y+=1.){
		double phase = 2.*M_PI*((double)rand()/(double)RAND_MAX);
		double xe =  amplitude_start/2.*cos(phase);
		double ye = -amplitude_start*sin(phase);
		particles[i].x  = xe;
		particles[i].y  = ye+y; 
		particles[i].z  = 0.;
		particles[i].vx = -2.*xe*OMEGA;
		particles[i].vy = 0.5*ye*OMEGA;
		particles[i].vz = 0;
		particles[i].ax = 0;
		particles[i].ay = 0;
		particles[i].az = 0;
		particles[i].m = g*g*G*OMEGA*OMEGA;
		i++;
	}
	// Do use ghost boxes in x and y
	nghostx = 0;
	nghosty = 4;
	nghostz = 0;
	cutoff_radius = ((double)nghosty-1.)*boxsize*((double)N-.5)/(double)N;
}

void problem_inloop(){
	int i=0;
	for (double y =-boxsize/2.0+boxsize/(double)N*0.5;i<N;y+=1.){
		struct particle p = particles[i];
		double e2 = (p.vx*p.vx/OMEGA/OMEGA+(2./OMEGA*p.vy+3.*p.x)*(2./OMEGA*p.vy+3.*p.x));
		if (e2>amplitude_finish*amplitude_finish){
			tmax =t;
		}
		/*
		if (p.x>.1||p.x<-.1){
			// Instability
			tmax = t;
		}
		if (p.y-y>.1||p.y-y<-.1){
			// Instability
			tmax = t;
		}
		*/
		i++;
	}
}

void problem_output(){

}

void problem_finish(){
	FILE *ofp;
	ofp = fopen("killtime.txt","a+"); 
	fprintf(ofp,"%f\t%f\n",g,t);
	fclose(ofp);
}
