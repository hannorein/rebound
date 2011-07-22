#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "tree.h"
#include "communication_mpi.h"

struct particle* 	particles;

int N_memory 		= 0;
int N 			= 0;
int Nmax		= 0;
int N_active_last 	= -1; // -1 is equivalent to N
int N_active_first 	= 0;

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
		fprintf(stderr,"ERROR: Number of root boxes (%d) not a multiple of mpi nodes (%d).\n",root_n,mpi_num);
		exit(-1);
	}
#endif // MPI

	boxsize_max = boxsize_x;
	if (boxsize_max<boxsize_y) boxsize_max = boxsize_y;
	if (boxsize_max<boxsize_z) boxsize_max = boxsize_z;
	
	printf("Initialized %d*%d*%d boxes.\n",root_nx,root_ny,root_nz);
}

void add_particle_local(struct particle pt);
void add_particle(struct particle pt){
#ifdef MPI
	int rootbox = get_rootbox_for_particle(pt);
	int root_n_per_node = root_n/mpi_num;
	int proc_id = rootbox/root_n_per_node;
	if (proc_id != mpi_id){
		// Add particle to array and send them to proc_id later. 
		add_particle_remote(pt,proc_id);
		return;
	}
#endif // MPI
	// Add particle to local partical array.
	add_particle_local(pt);
}


void add_particle_local(struct particle pt){
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


int get_rootbox_for_particle(struct particle pt){
	int i = ((int)floor((pt.x + boxsize_x/2.)/boxsize)+root_nx)%root_nx;
	int j = ((int)floor((pt.y + boxsize_y/2.)/boxsize)+root_ny)%root_ny;
	int k = ((int)floor((pt.z + boxsize_z/2.)/boxsize)+root_nz)%root_nz;
	int index = (k*root_ny+j)*root_nx+i;
	return index;
}

int get_rootbox_for_particle_int(int pt){
	struct particle p = particles[pt];
	int i = ((int)floor((p.x + boxsize_x/2.)/boxsize)+root_nx)%root_nx;
	int j = ((int)floor((p.y + boxsize_y/2.)/boxsize)+root_ny)%root_ny;
	int k = ((int)floor((p.z + boxsize_z/2.)/boxsize)+root_nz)%root_nz;
	int index = (k*root_ny+j)*root_nx+i;
	return index;
}
