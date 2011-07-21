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
#ifdef OPENGL
#include "opengl.h"
#endif // OPENGL
#ifdef MPI
#include "mpi.h"
void mpi_init(int argc, char** argv);
#endif // MPI
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
		printf("\nComputation finished. Total runtime: %f s\n",timing_final-timing_initial);
		exit(0);
	}

}

int main(int argc, char* argv[]) {
#ifdef MPI
	mpi_init(argc,argv);
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
	init_box();
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

#ifdef MPI
MPI_Datatype mpi_particle;
#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
MPI_Datatype mpi_cell;
#endif
MPI_Status stat; 
int mpi_num;
int mpi_id;

void mpi_init(int argc, char** argv){
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_num);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_id);
	
	
	// Setup MPI description of the particle structure 
	struct particle p;
	int bnum = 0;
	int blen[4];
	MPI_Aint indices[4];
	MPI_Datatype oldtypes[4];
#ifdef COLLISIONS_NONE
	blen[bnum] 	= 13;
#else //COLLISIONS_NONE
	blen[bnum] 	= 15; 
#endif //COLLISIONS_NONE
	indices[bnum] 	= 0; 
	oldtypes[bnum] 	= MPI_DOUBLE;
	bnum++;
#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
	blen[bnum] 	= 1; 
	indices[bnum] 	= (void*)&p.c - (void*)&p; 
	oldtypes[bnum] 	= MPI_CHAR;
	bnum++;
#endif
	blen[bnum] 	= 1; 
	indices[bnum] 	= sizeof(struct particle); 
	oldtypes[bnum] 	= MPI_UB;
	bnum++;
	MPI_Type_struct(bnum, blen, indices, oldtypes, &mpi_particle );

#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
	// Setup MPI description of the cell structure 
	struct cell c;
	bnum = 0;
#ifdef GRAVITY_TREE
#ifdef QUADRUPOLE
	blen[bnum] 	= 14;
#else // QUADRUPOLE
	blen[bnum] 	= 8;
#endif // QUADRUPOLE
#else //GRAVITY_TREE
	blen[bnum] 	= 4; 
#endif //GRAVITY_TREE
	indices[bnum] 	= 0; 
	oldtypes[bnum] 	= MPI_DOUBLE;
	bnum++;
	blen[bnum] 	= 1; 
	indices[bnum] 	= (void*)&c.oct - (void*)&c; 
	oldtypes[bnum] 	= MPI_CHAR;
	bnum++;
	blen[bnum] 	= 1; 
	indices[bnum] 	= (void*)&c.pt - (void*)&c; 
	oldtypes[bnum] 	= MPI_INT;
	bnum++;
	blen[bnum] 	= 1; 
	indices[bnum] 	= sizeof(struct cell); 
	oldtypes[bnum] 	= MPI_UB;
	bnum++;
	MPI_Type_struct(bnum, blen, indices, oldtypes, &mpi_cell );
#endif 

	if (mpi_id==0){
		printf("Using MPI with %d processors.\n",mpi_num);
	}

}
#endif // MPI
