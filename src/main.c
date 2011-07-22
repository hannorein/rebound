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
#include "tree.h"
#include "communication_mpi.h"
#ifdef OPENGL
#include "opengl.h"
#endif // OPENGL
#ifdef OPENMP
#include <omp.h>
#endif

double softening = 0.01;
double G=1;
double t=0;
double tmax=0; // Run forever
double dt = 0.001;

double timing_initial = -1;

void iterate(){	
	integrate_particles();
	t+=dt;  // Note: This might be better at the end of this function (t is used in check_boundaries).
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
		printf("\nComputation finished. Total runtime: %f s\n",timing_final-timing_initial);
		exit(0);
	}

}

int main(int argc, char* argv[]) {
#ifdef MPI
	communication_mpi_init(argc,argv);
#endif
#ifdef OPENMP
	printf("Using OpenMP with %d threads per node.\n",omp_get_num_threads());
#endif
	// Timing
	struct timeval tim;
	gettimeofday(&tim, NULL);
	timing_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
	// Initialiase random numbers, problem, box and OpengL
	srand ( time(NULL) );
	problem_init(argc, argv);
	boundaries_check();
	problem_output();
#ifdef OPENGL
	init_display(argc, argv);
#else // OPENGL
	while(1){
		iterate();
	}
#endif // OPENGL
	return 0;
}

