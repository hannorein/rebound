/**
 * @file 	main.c
 * @brief 	Main routine, iteration loop, timing.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <sys/time.h>
#include "integrator.h"
#include "boundaries.h"
#include "gravity.h"
#include "problem.h"
#include "output.h"
#include "collisions.h"
#include "tree.h"
#include "particle.h"
#include "main.h"
#include "communication_mpi.h"
#ifdef OPENGL
#include "display.h"
#endif // OPENGL
#ifdef OPENMP
#include <omp.h>
#endif
#ifdef GRAVITY_GRAPE
void gravity_finish();
#endif // GRAVITY_GRAPE

double softening 	= 0;
double G		= 1;
double t		= 0;
double tmax		= 0; 
double dt 		= 0.001;
double timing_initial 	= -1;
int    exit_simulation	= 0;

double boxsize 		= -1;
double boxsize_x 	= -1;
double boxsize_y 	= -1;
double boxsize_z 	= -1;
double boxsize_max 	= -1;
int root_nx		= 1;
int root_ny		= 1;
int root_nz		= 1;
int root_n		= 1;

void (*problem_additional_forces) () = NULL;

static char* 	logo[];		/**< Logo of rebound. */

void init_boxwidth(double _boxwidth){
	boxsize = _boxwidth;
	init_box();
}
void init_box(){	
	if (boxsize<=0 ){
		fprintf(stderr,"ERROR: Size of boxsize has to be set and be positive.\n");
		exit(-1);
	}
	if (root_nx <=0 || root_ny <=0 || root_nz <= 0){
		fprintf(stderr,"ERROR: Number of root boxes must be greater or equal to 1 in each direction.\n");
		exit(-1);
	}

	// Remove all particles
	free(particles);
	particles = NULL;

	// Setup box sizes
	boxsize_x = boxsize *(double)root_nx;
	boxsize_y = boxsize *(double)root_ny;
	boxsize_z = boxsize *(double)root_nz;
	root_n = root_nx*root_ny*root_nz;

#ifdef MPI
	// Make sure domain can be decomposed into equal number of root boxes per node.
	if ((root_n/mpi_num)*mpi_num != root_n){
		if (mpi_id==0) fprintf(stderr,"ERROR: Number of root boxes (%d) not a multiple of mpi nodes (%d).\n",root_n,mpi_num);
		exit(-1);
	}
#endif // MPI

	boxsize_max = boxsize_x;
	if (boxsize_max<boxsize_y) boxsize_max = boxsize_y;
	if (boxsize_max<boxsize_z) boxsize_max = boxsize_z;
	
#ifdef MPI
	printf("Initialized %d*%d*%d root boxes. MPI-node: %d. Process id: %d.\n",root_nx,root_ny,root_nz,mpi_id, getpid());
#else // MPI
	printf("Initialized %d*%d*%d root boxes. Process id: %d.\n",root_nx,root_ny,root_nz, getpid());
#endif // MPI
}

void iterate(){	
	// A 'DKD'-like integrator will do the first 'D' part.
	PROFILING_START()
	integrator_part1();
	PROFILING_STOP(PROFILING_CAT_INTEGRATOR)

	// Check for root crossings.
	PROFILING_START()
	boundaries_check();     
	PROFILING_STOP(PROFILING_CAT_BOUNDARY)

	// Update and simplify tree. 
	// Prepare particles for distribution to other nodes. 
	// This function also creates the tree if called for the first time.
	PROFILING_START()
#ifdef TREE
	tree_update();          
#endif //TREE

#ifdef MPI
	// Distribute particles and add newly received particles to tree.
	communication_mpi_distribute_particles();
#endif // MPI

#ifdef GRAVITY_TREE
	// Update center of mass and quadrupole moments in tree in preparation of force calculation.
	tree_update_gravity_data(); 
	
#ifdef MPI
	// Prepare essential tree (and particles close to the boundary needed for collisions) for distribution to other nodes.
	tree_prepare_essential_tree_for_gravity();

	// Transfer essential tree and particles needed for collisions.
	communication_mpi_distribute_essential_tree_for_gravity();
#endif // MPI
#endif // GRAVITY_TREE

	// Calculate accelerations. 
	gravity_calculate_acceleration();
	// Calculate non-gravity accelerations. 
	if (problem_additional_forces) problem_additional_forces();
	PROFILING_STOP(PROFILING_CAT_GRAVITY)

	// Call problem specific function. 
	problem_inloop();

	// A 'DKD'-like integrator will do the 'KD' part.
	PROFILING_START()
	integrator_part2();
	PROFILING_STOP(PROFILING_CAT_INTEGRATOR)

	// Do collisions here. We need both the positions and velocities at the same time.
#ifndef COLLISIONS_NONE
	// Check for root crossings.
	PROFILING_START()
	boundaries_check();     
	PROFILING_STOP(PROFILING_CAT_BOUNDARY)

	// Search for collisions using local and essential tree.
	PROFILING_START()
	collisions_search();

	// Resolve collisions (only local particles are affected).
	collisions_resolve();
	PROFILING_STOP(PROFILING_CAT_COLLISION)
#endif  // COLLISIONS_NONE

#ifdef OPENGL
	PROFILING_START()
	display();
	PROFILING_STOP(PROFILING_CAT_VISUALIZATION)
#endif // OPENGL
	problem_output();
	// Check if the simulation finished.
#ifdef MPI
	int _exit_simulation = 0;
	MPI_Allreduce(&exit_simulation, &_exit_simulation,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
	exit_simulation = _exit_simulation;
#endif // MPI
	// @TODO: Adjust timestep so that t==tmax exaclty at the end.
	if((t+dt>tmax && tmax!=0.0) || exit_simulation==1){
#ifdef GRAVITY_GRAPE
		gravity_finish();
#endif // GRAVITY_GRAPE
		problem_finish();
		struct timeval tim;
		gettimeofday(&tim, NULL);
		double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
		printf("\nComputation finished. Total runtime: %f s\n",timing_final-timing_initial);
#ifdef MPI
		MPI_Finalize();
#endif // MPI
		exit(0);
	}
}

int interrupt_counter = 0;
void interruptHandler(int var) {
	// This will try to quit the simulation nicely
	// at the end of the current timestep.
	switch(interrupt_counter){
		case 0:
			printf("\nInterrupt received. Will try to exit.\n");
			exit_simulation=1;
			break;
		default:
			printf("\nInterrupt received. Will exit immediately.\n");
			exit(-1);
	}
	interrupt_counter++;
}

int main(int argc, char* argv[]) {
#ifdef MPI
	communication_mpi_init(argc,argv);
	// Print logo only on main node.
	if (mpi_id==0){
#endif // MPI
		int i=0;
		while (logo[i]!=NULL){ printf("%s",logo[i++]); }
#ifdef MPI
		printf("Using MPI with %d nodes.\n",mpi_num);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif // MPI
#ifdef OPENMP
	printf("Using OpenMP with %d threads per node.\n",omp_get_max_threads());
#endif // OPENMP
	// Store time to calculate total runtime.
	struct timeval tim;
	gettimeofday(&tim, NULL);
	timing_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
	// Initialiase interrupts, random numbers, problem, box and OpengL
	signal(SIGINT, interruptHandler);
	signal(SIGKILL, interruptHandler);
	srand ( tim.tv_usec + getpid());
	problem_init(argc, argv);
	problem_output();
#ifdef OPENGL
	display_init(argc, argv);
#else // OPENGL
	while(1){
		// Main run loop.
		iterate();
	}
#endif // OPENGL
	return 0;
}


static char* logo[] = {
"                _                           _  \n",   
"               | |                         | | \n",  
"       _ __ ___| |__   ___  _   _ _ __   __| | \n", 
"      | '__/ _ \\ '_ \\ / _ \\| | | | '_ \\ / _` | \n", 
"      | | |  __/ |_) | (_) | |_| | | | | (_| | \n", 
"      |_|  \\___|_.__/ \\___/ \\__,_|_| |_|\\__,_| \n", 
"                                               \n",   
"                    `-:://::.`                 \n",
"                `/oshhoo+++oossso+:`           \n", 
"             `/ssooys++++++ossssssyyo:`        \n", 
"           `+do++oho+++osssso++++++++sy/`      \n", 
"          :yoh+++ho++oys+++++++++++++++ss.     \n", 
"         /y++hooyyooshooo+++++++++++++++oh-    \n", 
"        -dsssdssdsssdssssssssssooo+++++++oh`   \n", 
"        ho++ys+oy+++ho++++++++oosssssooo++so   \n", 
"       .d++oy++ys+++oh+++++++++++++++oosssod   \n", 
"       -h+oh+++yo++++oyo+++++++++++++++++oom   \n", 
"       `d+ho+++ys+++++oys++++++++++++++++++d   \n", 
"        yys++++oy+++++++oys+++++++++++++++s+   \n", 
"        .m++++++h+++++++++oys++++++++++++oy`   \n", 
"         -yo++++ss++++++++++oyso++++++++oy.    \n", 
"          .ss++++ho+++++++++++osys+++++yo`     \n", 
"            :ss+++ho+++++++++++++osssss-       \n", 
"              -ossoys++++++++++++osso.         \n", 
"                `-/oyyyssosssyso+/.            \n", 
"                      ``....`                  \n", 
"                                               \n",   
"    Copyright (C) 2011 Hanno Rein, Shangfei Liu\n",  
"       http://github.com/hannorein/rebound/    \n",    
"       http://arxiv.org/abs/1110.4876          \n",    
"                                               \n", 
NULL};
