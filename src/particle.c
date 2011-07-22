#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "tree.h"
#include "communication_mpi.h"

struct particle* 	particles;	// Main data stucture. Contains all particles.

int N 			= 0;		// Current number of particles on this node.
int Nmax		= 0;		// Current maximum number of particles in particle array.
int N_active_last 	= -1; 		// Last active (massive) particle. -1 is equivalent to N.
int N_active_first 	= 0;		// First active particle.


void particles_add_local(struct particle pt){
	if (Nmax<=N){
		Nmax += 128;
		particles = realloc(particles,sizeof(struct particle)*Nmax);
	}
	particles[N] = pt;
#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
	tree_add_particle_to_tree(N);
#endif
	N++;
}

void particles_add(struct particle pt){
#ifdef MPI
	int rootbox = particles_get_rootbox_for_particle(pt);
	int root_n_per_node = root_n/mpi_num;
	int proc_id = rootbox/root_n_per_node;
	if (proc_id != mpi_id){
		// Add particle to array and send them to proc_id later. 
		communication_mpi_add_particle_to_send_queue(pt,proc_id);
		return;
	}
#endif // MPI
	// Add particle to local partical array.
	particles_add_local(pt);
}




int particles_get_rootbox_for_particle(struct particle pt){
	int i = ((int)floor((pt.x + boxsize_x/2.)/boxsize)+root_nx)%root_nx;
	int j = ((int)floor((pt.y + boxsize_y/2.)/boxsize)+root_ny)%root_ny;
	int k = ((int)floor((pt.z + boxsize_z/2.)/boxsize)+root_nz)%root_nz;
	int index = (k*root_ny+j)*root_nx+i;
	return index;
}

