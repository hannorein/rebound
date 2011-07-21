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
	double offset = -0.5*boxsize_y + fmod(t*1.5*OMEGA*boxsize_x,boxsize_y);
#pragma omp parallel for
	for (int i=0;i<N;i++){
		// Radial
		while(particles[i].x>boxsize_x/2.){
			particles[i].x -= boxsize_x;
			particles[i].y += offset;
			particles[i].vy += 3./2.*OMEGA*boxsize_x;
		}
		while(particles[i].x<-boxsize_x/2.){
			particles[i].x += boxsize_x;
			particles[i].y -= offset;
			particles[i].vy -= 3./2.*OMEGA*boxsize_x;
		}
		// Azimuthal
		while(particles[i].y>boxsize_y/2.){
			particles[i].y -= boxsize_y;
		}
		while(particles[i].y<-boxsize_y/2.){
			particles[i].y += boxsize_y;
		}
		// Vertical (there should be no boundary, but periodic makes life easier)
		while(particles[i].z>boxsize_z/2.){
			particles[i].z -= boxsize_z;
		}
		while(particles[i].z<-boxsize_z/2.){
			particles[i].z += boxsize_z;
		}
	}
}

struct ghostbox get_ghostbox(int i, int j, int k){
	struct ghostbox gb;
	gb.shiftvx = 0;
	gb.shiftvy = -1.5*(double)i*OMEGA*boxsize_x;
	gb.shiftvz = 0;
	double shift = -fmod(gb.shiftvy*t,boxsize_y); 
	gb.shiftx = boxsize_x*(double)i;
	gb.shifty = boxsize_y*(double)j-shift;
	gb.shiftz = boxsize_z*(double)k;
	if(i>0) gb.shifty += boxsize_y*0.5;
	if(i<0) gb.shifty -= boxsize_y*0.5;
	return gb;
}


