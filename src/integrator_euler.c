#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"

void integrate_particles(){
	for (int i=0;i<N;i++){
		particles[i].vx += dt * particles[i].ax;
		particles[i].vy += dt * particles[i].ay;
		particles[i].vz += dt * particles[i].az;
		particles[i].x  += dt * particles[i].vx;
		particles[i].y  += dt * particles[i].vy;
		particles[i].z  += dt * particles[i].vz;
	}
}
	

