#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "integrator.h"
#include "main.h"
#include "boundaries.h"

int nghostx = 0;
int nghosty = 0;
int nghostz = 0;

void check_boundaries(){
	for (int i=0;i<N;i++){
		int removep = 0;
		if(particles[i].x>boxsize/2.){
			removep = 1;
		}
		if(particles[i].x<-boxsize/2.){
			removep = 1;
		}
		if(particles[i].y>boxsize/2.){
			removep = 1;
		}
		if(particles[i].y<-boxsize/2.){
			removep = 1;
		}
		if(particles[i].z>boxsize/2.){
			removep = 1;
		}
		if(particles[i].z<-boxsize/2.){
			removep = 1;
		}
		if (removep==1){
			// Note this has to be modified to work with the tree code.
			if (N==1){
				printf("Last particle removed. Exiting.");
				exit(0);
			}
			particles[i] = particles[N-1];
			i--;
			N--;
		}
	}
}

struct ghostbox get_ghostbox(int i, int j, int k){
	struct ghostbox gb;
	gb.shiftx = boxsize*(double)i;
	gb.shifty = boxsize*(double)j;
	gb.shiftz = boxsize*(double)k;
	gb.shiftvx = 0;
	gb.shiftvy = 0;
	gb.shiftvz = 0;
	return gb;
}


