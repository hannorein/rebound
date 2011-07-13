#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "opengl.h"

extern double OMEGA;
extern double coefficient_of_restitution;
extern double minimum_collision_velocity;
double cutoff_radius;
double g = 1;
double rp = 1;
double amplitude_start = 1e-5;
double amplitude_finish = 1e-5;

void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA = 1.;
	boxsize = 5; // = Number of particles
	dt = 1.12342535151342e-4*2.*M_PI;
	softening = 0;
	coefficient_of_restitution = 0.95;
	minimum_collision_velocity = .5;
	tmax = 100.*2.*M_PI;
	if (argc>=3){
		g  = atof(argv[1]);
		rp = atof(argv[2]);
	}else{
		fprintf(stderr,"No g and rp given!\n");
		exit(-1);
	}
	printf("Roche radius: %f\n",pow(g*g/3.,1./3.));
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
		particles[i].r = rp;
		i++;
	}
	// Do use ghost boxes in x and y
	nghostx = 0;
	nghosty = 2;
	nghostz = 0;
	cutoff_radius = ((double)nghosty-1.)*boxsize*((double)N-.5)/(double)N;
}

void problem_inloop(){
}

void output_mean_abs_position(char* filename){
	// Algorithm with reduced roundoff errors (see wikipedia)
	double x=0,y=0,z=0;
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		x = x + (fabs(p.x) - x)/(double)(i+1);
		y = y + (fabs(p.y) - y)/(double)(i+1);
		z = z + (fabs(p.z) - z)/(double)(i+1);
	}
	FILE* of = fopen(filename,"a"); 
	fprintf(of,"%e\t%e\t%e\t%e\n",t,x,y,z);
	fclose(of);
}

void output_mean_epicycle_amplitude(char* filename){
	// Algorithm with reduced roundoff errors (see wikipedia)
	double e=0;
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		double ep = sqrt(p.vx*p.vx/(OMEGA*OMEGA)+(2./OMEGA*p.vy+3.*p.x)*(2./OMEGA*p.vy+3.*p.x));
		e = e + (ep - e)/(double)(i+1);
	}
	FILE* of = fopen(filename,"a"); 
	fprintf(of,"%e\t%e\n",t,e);
	fclose(of);
}

int pngoutputnum = 0;
void problem_output(){
	if (output_check(1e-2*2.*M_PI)){
		//char buf[1024];
		//sprintf(buf,"veldisp_g%0.3f_rp%0.3f.txt",g,rp);
		//output_mean_abs_position(buf);
		output_mean_epicycle_amplitude("mean_epicycle_amplitude.txt");
	}
#ifdef OPENGL
#ifdef LIBPNG
	if (output_check(2e-2*2.*M_PI)&&display_init_done){
		char buf[1024];
		sprintf(buf,"png/%09d.png",pngoutputnum);
		output_png(buf);
		pngoutputnum++;
	}
#endif
#endif
}

void problem_finish(){
}
