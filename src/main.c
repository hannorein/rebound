#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "integrator.h"
#include "boundaries.h"
#include "gravity.h"
#ifdef OPENGL
#include "opengl.h"
#endif

double boxsize = 1;
double softening = 0.01;
double G=1;
double t=0;
double tmax=10;
double dt = 0.01;
int N = 1000;

void iterate(){	
	integrate_particles();
	t+=dt;
	check_boundaries();
	printf("t = %f\n",t);
#ifdef OPENGL
	display();
#endif
}

int main(int argc, char* argv[]) {
	srand ( time(NULL) );
	init_particles();
#ifdef OPENGL
	init_display(argc, argv);
#else
	while(t<tmax){
		iterate();
	}
#endif
	printf("Computation finished\n");
	return 0;
}
	
