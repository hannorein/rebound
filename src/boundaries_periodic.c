#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "integrator.h"
#include "main.h"

void check_boundaries(){
	for (int i=0;i<N;i++){
		while(particles[i].x>boxsize/2.){
			particles[i].x -= boxsize;
		}
		while(particles[i].x<-boxsize/2.){
			particles[i].x += boxsize;
		}
		while(particles[i].y>boxsize/2.){
			particles[i].y -= boxsize;
		}
		while(particles[i].y<-boxsize/2.){
			particles[i].y += boxsize;
		}
		while(particles[i].z>boxsize/2.){
			particles[i].z -= boxsize;
		}
		while(particles[i].z<-boxsize/2.){
			particles[i].z += boxsize;
		}
	}
}
