#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"

struct particle* restrict particles;

void init_particles(int _N){	
	N = _N;
	N_active_first = 0;
	N_active_last  = _N;
	if (boxsize!=-1){
		boxsize_x = boxsize;
		boxsize_y = boxsize;
		boxsize_z = boxsize;
	}
	particles = calloc(N,sizeof(struct particle));
	printf("Initialized memory for %d particles.\n",N);
}
