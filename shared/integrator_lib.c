#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "particle.h"
#include "integrator.h"
#include "gravity.h"
#include "particle.h"
#include "main.h"

double dt 	= 0.01;	// Default values
double t 	= 0;
double G 	= 1;
double softening = 0;

void (*problem_additional_forces) () = NULL;

struct particle* particles = NULL;
void setp(struct particle* _p){
	free(particles);
	particles = malloc(sizeof(struct particle)*N);
	memcpy(particles,_p,sizeof(struct particle)*N);
}
struct particle particle_get(int i){
	return particles[i];
}

void step(){
	if (N<=0){
		fprintf(stderr,"\n\033[1mError!\033[0m No particles found. Exiting.\n");
		return;
	}
	integrator_part1();
	gravity_calculate_acceleration();
	integrator_part2();
}

void integrate(double tmax){
	while(t<tmax){
		if (N<=0){
			fprintf(stderr,"\n\033[1mError!\033[0m No particles found. Exiting.\n");
			return;
		}
		step();
		if (t+dt>=tmax){
			dt = tmax-t;
		}
	}
}
	 

// Infinite box size
int root_nx = 1;
int root_ny = 1;
int root_nz = 1;
double boxsize = 0;
double boxsize_x = 0;
double boxsize_y = 0;
double boxsize_z = 0;

int boundaries_particle_is_in_box(struct particle p){
	return 1;
}

// No ghost boxes for now.
int nghostx = 0;	
int nghosty = 0;
int nghostz = 0;
struct ghostbox{
	double shiftx;		/**< Relative x position */
	double shifty;		/**< Relative y position */
	double shiftz;		/**< Relative z position */
	double shiftvx;		/**< Relative x velocity */
	double shiftvy;		/**< Relative y velocity */
	double shiftvz;		/**< Relative z velocity */
};

struct ghostbox boundaries_get_ghostbox(int i, int j, int k){
	struct ghostbox gb;
	gb.shiftx = 0;
	gb.shifty = 0;
	gb.shiftz = 0;
	gb.shiftvx = 0;
	gb.shiftvy = 0;
	gb.shiftvz = 0;
	return gb;
}
