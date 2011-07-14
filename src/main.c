#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "integrator.h"
#include "boundaries.h"
#include "gravity.h"
#include "problem.h"
#include "collisions.h"
#ifdef OPENGL
#include "opengl.h"
#endif

double boxsize = 1;
double softening = 0.01;
double G=1;
double t=0;
double tmax=0; // Run forever
double dt = 0.001;
int N = 0;
int N_active = 0;

void iterate(){	
	integrate_particles();
	t+=dt;
	check_boundaries();
	printf("t = %f\n",t);
#ifdef OPENGL
	display();
#endif
	collisions_search();
	collisions_resolve();
	problem_inloop();
	problem_output();
	while(t>=tmax && tmax!=0.0){
		problem_finish();
		exit(0);
	}

}

int main(int argc, char* argv[]) {
	srand ( time(NULL) );
	problem_init(argc, argv);
	problem_output();
#ifdef OPENGL
	init_display(argc, argv);
#else
	while(t<tmax||tmax==0.0){
		iterate();
	}
#endif
	printf("Computation finished\n");
	return 0;
}
	
