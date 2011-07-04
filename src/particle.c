#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"

struct particle* particles;

void init_particles(){
	printf("Initializing %d particles\n",N);
	particles = malloc(sizeof(struct particle)*N);
	for (int i =0;i<N;i++){
		particles[i].x = ((double)rand()/(double)RAND_MAX-0.5)*boxsize;
		particles[i].y = ((double)rand()/(double)RAND_MAX-0.5)*boxsize;
		particles[i].z = ((double)rand()/(double)RAND_MAX-0.5)*boxsize;
		particles[i].vx = 0;
		particles[i].vy = 0;
		particles[i].vz = 0;
		particles[i].ax = 0;
		particles[i].ay = 0;
		particles[i].az = 0;
		particles[i].m = 0.001;
	}
}
	 
