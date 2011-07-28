#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "integrator.h"
#include "main.h"
#include "boundaries.h"
#include "tree.h"

int nghostx = 0;
int nghosty = 0;
int nghostz = 0;

void boundaries_check(){
	for (int i=0;i<N;i++){
		int removep = 0;
		if(particles[i].x>boxsize_x/2.){
			removep = 1;
		}
		if(particles[i].x<-boxsize_x/2.){
			removep = 1;
		}
		if(particles[i].y>boxsize_y/2.){
			removep = 1;
		}
		if(particles[i].y<-boxsize_y/2.){
			removep = 1;
		}
		if(particles[i].z>boxsize_z/2.){
			removep = 1;
		}
		if(particles[i].z<-boxsize_z/2.){
			removep = 1;
		}
		if (removep==1){
			if (N==1){
				printf("Last particle removed. Exiting.\n");
				exit(0);
			}
			// Note this has to be modified to work with the tree code.
#ifdef TREE
#error OPEN_BOUNDARY CONDITIONS NOT YET SUPPORTYED WITH TREE
			printf("Particle leaving the box. Exiting.\n");
			exit(0);
#else
			particles[i] = particles[N-1];
			i--;
			N--;
#endif
		}
	}
}

struct ghostbox boundaries_get_ghostbox(int i, int j, int k){
	struct ghostbox gb;
	gb.shiftx = boxsize_x*(double)i;
	gb.shifty = boxsize_y*(double)j;
	gb.shiftz = boxsize_z*(double)k;
	gb.shiftvx = 0;
	gb.shiftvy = 0;
	gb.shiftvz = 0;
	return gb;
}


