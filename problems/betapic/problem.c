#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"

void read_bb();
double get_bb(double a);
double* bb_a=NULL;
double* bb=NULL;
int bb_N=0;

void problem_init(int argc, char* argv[]){
	// Setup constants
	boxsize_x = 5; 
	boxsize_y = 5; 
	boxsize_z = 5; 
	dt = 2.0e-3*2.*M_PI;
	// Setup particle structures
	init_particles(10000); // Number of particles
	N_active_first = 1; // Only the planet's gravity is felt by all other particles
	N_active_last = 2; // Only the planet's gravity is felt by all other particles
	// Initial conditions
	// Planet
	double planet_e = 0.2;
	particles[0].x  = 0;
	particles[0].y  = 0; 
	particles[0].vx = 0;
	particles[0].vy = 0;
	particles[0].m  = 1;
	particles[1].x  = 1.-planet_e;
	particles[1].y  = 0; 
	particles[1].vx = 0;
	particles[1].vy = sqrt(2./(1.-planet_e)-1.);
	particles[1].m  = 4.57e-3;
	// Test particles
	read_bb();
	printf("Loaded %d values from bb file.\n",bb_N);


	int i=2;
	while(i<N){
		double x = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_x;
		double y = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
		double a = sqrt(x*x+y*y);
		x = a; y=0;
		//double phi = atan2(y,x);

		double e = planet_e*get_bb(a);
		if (a>boxsize_x/2.) continue;
		double vkep = sqrt((particles[0].m)/a);
		particles[i].x  = x*(1.-e);
		particles[i].y  = y; 
		particles[i].z  = 1e-2*vkep*((double)rand()/(double)RAND_MAX-0.5);
		particles[i].vx = 0;
		particles[i].vy = vkep*sqrt(2./(1.-e)-1.);
		i++;
	}
}



void read_bb(){
	FILE* inf = fopen("beta2beta1.tsv","r"); 
	double bb_t_a, bb_t;
	bb_N=0;
	while( fscanf(inf,"%lf\t%lf",&bb_t_a, &bb_t)!=EOF ){
		bb_a = realloc(bb_a, sizeof(double)*(bb_N+1));	
		bb = realloc(bb, sizeof(double)*(bb_N+1));	
		bb[bb_N] = bb_t;
		bb_a[bb_N] = bb_t_a;
		bb_N++;
	}
	fclose(inf);
	
}

double get_bb(double a){
	int i =1;
	while (bb_a[i]<a){
		if (i>=bb_N-2) break;
		i++;
	}
	double slope = (bb[i]-bb[i-1])/(bb_a[i]-bb_a[i-1]);
	return bb[i]+slope*( a - (bb_a[i]+bb_a[i-1])/2. ); 
}

int firstoutput=1;

void problem_output(){
	// Output every 1.23 orbits or so. Start after 50 orbits.
	if (output_check(1.23456*2.*M_PI)&&t>100.*M_PI){
		if (firstoutput==1){
			firstoutput=0;
			output_ascii("particles.txt");
		}else{
			output_append_ascii("particles.txt");
		}
	}
	if (output_check(2.*M_PI)){
		output_timing();
	}
}

void problem_finish(){
}

void problem_inloop(){
}
