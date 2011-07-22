#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "tree.h"
#ifdef MPI
#include "mpi.h"
#endif // MPI

struct particle* particles;
struct particle** particles_send;
int* particles_send_N;
int* particles_send_Nmax;
struct particle** particles_recv;
int* particles_recv_N;
int* particles_recv_Nmax;
void add_particle_remote(struct particle pt, int proc_id);
void add_particle_local(struct particle pt);

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
double boxsize_min 	= -1;

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
	
	particles_send   	= calloc(mpi_num,sizeof(struct particle*));
	particles_send_N 	= calloc(mpi_num,sizeof(int));
	particles_send_Nmax 	= calloc(mpi_num,sizeof(int));
	particles_recv   	= calloc(mpi_num,sizeof(struct particle*));
	particles_recv_N 	= calloc(mpi_num,sizeof(int));
	particles_recv_Nmax 	= calloc(mpi_num,sizeof(int));

#endif // MPI

	boxsize_max = boxsize_x;
	if (boxsize_max<boxsize_y) boxsize_max = boxsize_y;
	if (boxsize_max<boxsize_z) boxsize_max = boxsize_z;
	
	printf("Initialized %d*%d*%d boxes.\n",root_nx,root_ny,root_nz);
}

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

#ifdef MPI
void particles_communicate(){
	// Distribute the number of particles to be transferred.
	for (int i=0;i<mpi_num;i++){
		MPI_Scatter(particles_send_N, 1, MPI_INT, &(particles_recv_N[i]), 1, MPI_INT, i, MPI_COMM_WORLD);
	}
	// Allocate memory for incoming particles
	for (int i=0;i<mpi_num;i++){
		if  (i==mpi_id) continue;
		if (particles_recv_Nmax[i]<particles_recv_N[i]){
			particles_recv_Nmax[i] += 32;
			particles_recv[i] = realloc(particles_recv[i],sizeof(struct particle)*particles_recv_Nmax[i]);
		}
	}

	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (particles_recv_N[i]==0) continue;
	}
	
	// Exchange particles via MPI.
	// Using non-blocking receive call.
	MPI_Request request[mpi_num];
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (particles_recv_N[i]==0) continue;
		MPI_Irecv(particles_recv[i], particles_recv_N[i], mpi_particle, i, i*mpi_num+mpi_id, MPI_COMM_WORLD, &(request[i]));
	}
	// Using blocking send call.
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (particles_send_N[i]==0) continue;
		MPI_Send(particles_send[i], particles_send_N[i], mpi_particle, i, mpi_id*mpi_num+i, MPI_COMM_WORLD);
	}
	// Wait for all particles to be received.
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (particles_recv_N[i]==0) continue;
		MPI_Status status;
		MPI_Wait(&(request[i]), &status);
	}
	// Add particles to local tree
	for (int i=0;i<mpi_num;i++){
		for (int j=0;j<particles_recv_N[i];j++){
			add_particle(particles_recv[i][j]);
		}
	}
	// Bring everybody into sync, clean up. 
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<mpi_num;i++){
		particles_send_N[i] = 0;
		particles_recv_N[i] = 0;
	}
}


void add_particle_remote(struct particle pt, int proc_id){
	if (particles_send_Nmax[proc_id]<=particles_send_N[proc_id]){
		particles_send_Nmax[proc_id] += 128;
		particles_send[proc_id] = realloc(particles_send[proc_id],sizeof(struct particle)*particles_send_Nmax[proc_id]);
	}
	particles_send[proc_id][particles_send_N[proc_id]] = pt;
	particles_send_N[proc_id]++;
}
#endif // MPI

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
