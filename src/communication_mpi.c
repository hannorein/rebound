#ifdef MPI
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "particle.h"
#include "tree.h"
#include "communication_mpi.h"

#include "mpi.h"
MPI_Datatype mpi_particle;
#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
MPI_Datatype mpi_cell;
#endif
MPI_Status stat; 
int mpi_num;
int mpi_id;

struct particle** 	particles_send;
int* 			particles_send_N;
int* 			particles_send_Nmax;
struct particle** 	particles_recv;
int* 			particles_recv_N;
int* 			particles_recv_Nmax;

void communication_mpi_init(int argc, char** argv){
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_num);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_id);
	
	
	// Setup MPI description of the particle structure 
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
	struct particle p;
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
	MPI_Type_commit(&mpi_particle); 

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
	MPI_Type_commit(&mpi_cell); 
#endif 
	
	// Prepare send/recv buffers.
	particles_send   	= calloc(mpi_num,sizeof(struct particle*));
	particles_send_N 	= calloc(mpi_num,sizeof(int));
	particles_send_Nmax 	= calloc(mpi_num,sizeof(int));
	particles_recv   	= calloc(mpi_num,sizeof(struct particle*));
	particles_recv_N 	= calloc(mpi_num,sizeof(int));
	particles_recv_Nmax 	= calloc(mpi_num,sizeof(int));


	if (mpi_id==0){
		printf("Using MPI with %d nodes.\n",mpi_num);
	}

}

void communication_mpi_distribute_particles(){
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
			particles_add(particles_recv[i][j]);
		}
	}
	// Bring everybody into sync, clean up. 
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<mpi_num;i++){
		particles_send_N[i] = 0;
		particles_recv_N[i] = 0;
	}
}

void communication_mpi_add_particle_to_send_queue(struct particle pt, int proc_id){
	int send_N = particles_send_N[proc_id];
	if (particles_send_Nmax[proc_id] <= send_N){
		particles_send_Nmax[proc_id] += 128;
		particles_send[proc_id] = realloc(particles_send[proc_id],sizeof(struct particle)*particles_send_Nmax[proc_id]);
	}
	particles_send[proc_id][send_N] = pt;
	particles_send_N[proc_id]++;
}



#endif // MPI
