#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "boundaries.h"
#include "main.h"
#include "tree.h"

extern const double OMEGA;
int nghostx = 1;
int nghosty = 0;
int nghostz = 0;

void check_boundaries(){
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		// No boundary in the radial direction. This will not work with GRAVITY_TREE!
		
		// Azimuthal
		while(particles[i].y>boxsize/2.){
			particles[i].y -= boxsize;
		}
		while(particles[i].y<-boxsize/2.){
			particles[i].y += boxsize;
		}
		// Vertical (there should be no boundary, but periodic makes life easier)
		while(particles[i].z>boxsize/2.){
			particles[i].z -= boxsize;
		}
		while(particles[i].z<-boxsize/2.){
			particles[i].z += boxsize;
		}
	}
}

struct ghostbox boundaries_get_ghostbox(int i, int j, int k){
	struct ghostbox gb;
	gb.shiftvx = 0;
	gb.shiftvy = 0;
	gb.shiftvz = 0;
	gb.shiftx  = 0;
	gb.shifty  = boxsize*(double)j;
	gb.shiftz  = boxsize*(double)k;
	return gb;
}


