#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "integrator.h"
#include "main.h"
#include "boundaries.h"

int nghostx = 1;
int nghosty = 1;
int nghostz = 1;

void check_boundaries(){
	for (int i=0;i<N;i++){
		while(particles[i].x>boxsize_x/2.){
			particles[i].x -= boxsize_x;
		}
		while(particles[i].x<-boxsize_x/2.){
			particles[i].x += boxsize_x;
		}
		while(particles[i].y>boxsize_y/2.){
			particles[i].y -= boxsize_y;
		}
		while(particles[i].y<-boxsize_y/2.){
			particles[i].y += boxsize_y;
		}
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
	gb.shiftx = boxsize_x*(double)i;
	gb.shifty = boxsize_x*(double)j;
	gb.shiftz = boxsize_y*(double)k;
	gb.shiftvx = 0;
	gb.shiftvy = 0;
	gb.shiftvz = 0;
	return gb;
}


