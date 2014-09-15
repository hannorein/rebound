#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "particle.h"

int N 		= 0;
int N_active 	= -1; 	// All particles have mass
double dt 	= 0.1;	// Default values
double t 	= 0;
double G 	= 1;
double softening = 0;

void (*problem_additional_forces) () = NULL;

struct particle* particles = NULL;
void setp(struct particle* _p){
	particles = malloc(sizeof(struct particle)*N);
	memcpy(particles,_p,sizeof(struct particle)*N);
}
struct particle* getp(){
	return particles;
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
