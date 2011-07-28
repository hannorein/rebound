/**
 * @file 	main.c
 * @brief 	Main routine, iteration loop, timing.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of nbody.
 *
 * nbody is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * nbody is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with nbody.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
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
double tmax		= 0; 
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
	// Make sure domain can be decomposed into equal number of root boxes per node.
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
	// A 'DKD'-like integrator will do the first 'D' part.
	integrator_part1();

	// Check for root crossings.
	boundaries_check();     

	// Update and simplify tree. 
	// Prepare particles for distribution to other nodes. 
	// This function also creates the tree if called for the first time.
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

	// Call problem specific function (e.g. to add additional forces). 
	problem_inloop();

	// A 'DKD'-like integrator will do the 'KD' part.
	integrator_part2();

	// Do collisions here. We need both the positions and velocities at the same time.
#ifndef COLLISIONS_NONE
	// Check for root crossings.
	boundaries_check();     

#ifdef COLLISIONS_TREE
	// Update and simplify tree. 
	// Prepare particles for distribution to other nodes. 
	tree_update();          

#ifdef MPI
	// Distribute particles and add newly received particles to tree.
	communication_mpi_distribute_particles();
	
	// Prepare essential tree (and particles close to the boundary needed for collisions) for distribution to other nodes.
	tree_prepare_essential_tree_for_collisions();

	// Transfer essential tree and particles needed for collisions.
	communication_mpi_distribute_essential_tree_for_collisions();
#endif // MPI
#endif // COLLISIONS_TREE

	// Search for collisions using local and essential tree.
	collisions_search();

	// Resolve collisions (only local particles are affected).
	collisions_resolve();
#endif  // COLLISIONS_NONE

#ifdef OPENGL
	display();
#endif // OPENGL
	problem_output();
	// Check if the simulation finished.
	// @TODO: Adjust timestep so that t==tmax exaclty at the end.
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
	// Store time to calculate total runtime.
	struct timeval tim;
	gettimeofday(&tim, NULL);
	timing_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
	// Initialiase random numbers, problem, box and OpengL
	srand ( tim.tv_usec + getpid());
	problem_init(argc, argv);
	problem_output();
#ifdef OPENGL
	init_display(argc, argv);
#else // OPENGL
	while(1){
		// Main run loop.
		iterate();
	}
#endif // OPENGL
	return 0;
}

