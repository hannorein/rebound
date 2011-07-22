#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "gravity.h"

// Leapfrog integrator (Drift-Kick-Drift)
// for non-rotating frame.
void integrate_particles(){
	calculate_forces();
#pragma omp parallel for
	for (int i=0;i<N;i++){
		particles[i].x  += dt * particles[i].vx;
		particles[i].y  += dt * particles[i].vy;
		particles[i].z  += dt * particles[i].vz;
		particles[i].vx += dt * particles[i].ax;
		particles[i].vy += dt * particles[i].ay;
		particles[i].vz += dt * particles[i].az;
	}
	boundaries_check();
}
	

