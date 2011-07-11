#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "boundaries.h"
#include "main.h"

extern const double OMEGA;
int nghostx = 1;
int nghosty = 1;
int nghostz = 0;

void check_boundaries(){
	double offset = -0.5*boxsize + fmod(1.5*OMEGA*t,1.)*boxsize;
	for (int i=0;i<N;i++){
		// Radial
		while(particles[i].x>boxsize/2.){
			particles[i].x -= boxsize;
			particles[i].y += offset;
			particles[i].vy += 3./2.*OMEGA*boxsize;
		}
		while(particles[i].x<-boxsize/2.){
			particles[i].x += boxsize;
			particles[i].y -= offset;
			particles[i].vy -= 3./2.*OMEGA*boxsize;
		}
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

struct ghostbox get_ghostbox(int i, int j, int k){
	gb.shiftvx = 0;
	gb.shiftvy = 1.5*(double)i*OMEGA*boxsize;
	gb.shiftvz = 0;
	double shift = fmod(gb.shiftvy*t,1.)*boxsize; 
	struct ghostbox gb;
	gb.shiftx = boxsize*(double)i;
	gb.shifty = boxsize*(double)j-shift;
	gb.shiftz = boxsize*(double)k;
	if(i>0) gb.shifty += boxsize*0.5;
	if(i<0) gb.shifty -= boxsize*0.5;
	return gb;
}


