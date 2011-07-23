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

double softening 	= 0;
double G		= 1;
double t		= 0;
double tmax		= 0; // 0 == run forever
double dt 		= 0.001;
double timing_initial 	= -1;

double boxsize 		= -1;
double boxsize_x 	= -1;
double boxsize_y 	= -1;
double boxsize_z 	= -1;
double boxsize_max 	= -1;
int root_nx		= 1;
int root_ny		= 1;
int root_nz		= 1;
int root_n		= 1;

void init_box(){	
	boxsize_x = boxsize *(double)root_nx;
	boxsize_y = boxsize *(double)root_ny;
	boxsize_z = boxsize *(double)root_nz;
	root_n = root_nx*root_ny*root_nz;

#ifdef MPI
	if ((root_n/mpi_num)*mpi_num != root_n){
		if (mpi_id==0) fprintf(stderr,"ERROR: Number of root boxes (%d) not a multiple of mpi nodes (%d).\n",root_n,mpi_num);
		exit(-1);
	}
#endif // MPI

	boxsize_max = boxsize_x;
	if (boxsize_max<boxsize_y) boxsize_max = boxsize_y;
	if (boxsize_max<boxsize_z) boxsize_max = boxsize_z;
	
	printf("Initialized %d*%d*%d root boxes.\n",root_nx,root_ny,root_nz);
}


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
#endif // MPI
#ifdef OPENMP
	printf("Using OpenMP with %d threads per node.\n",omp_get_max_threads());
#endif // OPENMP
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

