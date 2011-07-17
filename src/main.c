#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "particle.h"
#include "integrator.h"
#include "boundaries.h"
#include "gravity.h"
#include "problem.h"
#include "collisions.h"
#ifdef OPENGL
#include "opengl.h"
#endif

double boxsize = -1;
double boxsize_x = -1;
double boxsize_y = -1;
double boxsize_z = -1;
double boxsize_max = -1;
double boxsize_min = -1;
double softening = 0.01;
double G=1;
double t=0;
double tmax=0; // Run forever
double dt = 0.001;
int N = 0;
int N_active_last = -1;
int N_active_first = -1;

double timing_initial = -1;

void iterate(){	
	integrate_particles();
	t+=dt;
	check_boundaries();
#ifdef OPENGL
	display();
#endif
	collisions_search();
	collisions_resolve();
	problem_inloop();
	problem_output();
	if(t+dt>tmax && tmax!=0.0){
		problem_finish();
		struct timeval tim;
		gettimeofday(&tim, NULL);
		double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
		printf("Computation finished. Total runtime: %f s\n",timing_final-timing_initial);
		exit(0);
	}

}

int main(int argc, char* argv[]) {
	// Timing
	struct timeval tim;
	gettimeofday(&tim, NULL);
	timing_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
	// Initialiase random numbers, problem, box and OpengL
	srand ( time(NULL) );
	problem_init(argc, argv);
	boxsize_max = boxsize_x;
	if (boxsize_max<boxsize_y) boxsize_max = boxsize_y;
	if (boxsize_max<boxsize_z) boxsize_max = boxsize_z;
	boxsize_min = boxsize_x;
	if (boxsize_y<boxsize_min) boxsize_min = boxsize_y;
	if (boxsize_z<boxsize_min) boxsize_min = boxsize_z;
	problem_output();
#ifdef OPENGL
	init_display(argc, argv);
#else
	while(1){
		iterate();
	}
#endif
	return 0;
}
	
