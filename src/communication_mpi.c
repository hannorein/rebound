#ifdef MPI
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "particle.h"
#include "main.h"
#include "tree.h"
#include "boundaries.h"
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

#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
struct cell** 	tree_essential_send;
int* 		tree_essential_send_N;
int* 		tree_essential_send_Nmax;
struct cell** 	tree_essential_recv;
int* 		tree_essential_recv_N;
int* 		tree_essential_recv_Nmax;
#endif


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
	
	// Prepare send/recv buffers for particles
	particles_send   	= calloc(mpi_num,sizeof(struct particle*));
	particles_send_N 	= calloc(mpi_num,sizeof(int));
	particles_send_Nmax 	= calloc(mpi_num,sizeof(int));
	particles_recv   	= calloc(mpi_num,sizeof(struct particle*));
	particles_recv_N 	= calloc(mpi_num,sizeof(int));
	particles_recv_Nmax 	= calloc(mpi_num,sizeof(int));

#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
	// Prepare send/recv buffers for essential tree
	tree_essential_send   	= calloc(mpi_num,sizeof(struct cell*));
	tree_essential_send_N 	= calloc(mpi_num,sizeof(int));
	tree_essential_send_Nmax= calloc(mpi_num,sizeof(int));
	tree_essential_recv   	= calloc(mpi_num,sizeof(struct cell*));
	tree_essential_recv_N 	= calloc(mpi_num,sizeof(int));
	tree_essential_recv_Nmax= calloc(mpi_num,sizeof(int));
#endif

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
		while (particles_recv_Nmax[i]<particles_recv_N[i]){
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
	while (particles_send_Nmax[proc_id] <= send_N){
		particles_send_Nmax[proc_id] += 128;
		particles_send[proc_id] = realloc(particles_send[proc_id],sizeof(struct particle)*particles_send_Nmax[proc_id]);
	}
	particles_send[proc_id][send_N] = pt;
	particles_send_N[proc_id]++;
}

#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)

struct aabb{ // axis aligned bounding box
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;
} aabb;

struct aabb communication_boundingbox_for_root(int index){
	int i = index%root_nx;
	int j = ((index-i)/root_nx)%root_ny;
	int k = ((index-i)/root_nx-j)/root_ny;
	struct aabb boundingbox;
	boundingbox.xmin = -boxsize_x/2.+boxsize*(double)i;
	boundingbox.ymin = -boxsize_y/2.+boxsize*(double)j;
	boundingbox.zmin = -boxsize_z/2.+boxsize*(double)k;
	boundingbox.xmax = -boxsize_x/2.+boxsize*(double)(i+1);
	boundingbox.ymax = -boxsize_y/2.+boxsize*(double)(j+1);
	boundingbox.zmax = -boxsize_z/2.+boxsize*(double)(k+1);
	return boundingbox;
}

struct aabb communication_boundingbox_for_proc(int proc_id){
	int root_n_per_node = root_n/mpi_num;
	int root_start = proc_id*root_n_per_node;
	int root_stop  = (proc_id+1)*root_n_per_node;
	struct aabb boundingbox = communication_boundingbox_for_root(root_start);
	for (int i=root_start+1;i<root_stop;i++){
		struct aabb boundingbox2 = communication_boundingbox_for_root(i);
		if (boundingbox.xmin > boundingbox2.xmin) boundingbox.xmin = boundingbox2.xmin;
		if (boundingbox.ymin > boundingbox2.ymin) boundingbox.ymin = boundingbox2.ymin;
		if (boundingbox.zmin > boundingbox2.zmin) boundingbox.zmin = boundingbox2.zmin;
		if (boundingbox.xmax < boundingbox2.xmax) boundingbox.xmax = boundingbox2.xmax;
		if (boundingbox.ymax < boundingbox2.ymax) boundingbox.ymax = boundingbox2.ymax;
		if (boundingbox.zmax < boundingbox2.zmax) boundingbox.zmax = boundingbox2.zmax;
	}
	return boundingbox;
}
double communication_distance2_of_aabb_to_cell(struct aabb bb, struct cell* node){
	double distancex = fabs(node->x - (bb.xmin+bb.xmax)/2.)  -  (node->w + bb.xmax-bb.xmin)/2.;
	double distancey = fabs(node->y - (bb.ymin+bb.ymax)/2.)  -  (node->w + bb.ymax-bb.ymin)/2.;
	double distancez = fabs(node->z - (bb.zmin+bb.zmax)/2.)  -  (node->w + bb.zmax-bb.zmin)/2.;
	if (distancex<=0 && distancey<=0 && distancez<=0) return 0;	// Overlapping
	if (distancex<=0 && distancey <=0) return distancez*distancez;
	if (distancey<=0 && distancez <=0) return distancex*distancex;
	if (distancez<=0 && distancex <=0) return distancey*distancey;
	if (distancex<=0)  return distancey*distancey + distancez*distancez;
	if (distancey<=0)  return distancex*distancex + distancez*distancez;
	if (distancez<=0)  return distancey*distancey + distancex*distancex;
	return distancex*distancex + distancey*distancey + distancez*distancez;
}

double communication_distance2_of_proc_to_node(int proc_id, struct cell* node){
	int nghostxcol = (nghostx>1?1:0);
	int nghostycol = (nghosty>1?1:0);
	int nghostzcol = (nghostz>1?1:0);
	double distance = boxsize*(double)root_n;
	for (int gbx=-nghostxcol; gbx<=nghostxcol; gbx++){
	for (int gby=-nghostycol; gby<=nghostycol; gby++){
	for (int gbz=-nghostzcol; gbz<=nghostzcol; gbz++){
		struct ghostbox gb = get_ghostbox(gbx,gby,gbz);
		struct aabb boundingbox = communication_boundingbox_for_proc(proc_id);
		boundingbox.xmin+=gb.shiftx;
		boundingbox.xmax+=gb.shiftx;
		boundingbox.ymin+=gb.shifty;
		boundingbox.ymax+=gb.shifty;
		boundingbox.zmin+=gb.shiftz;
		boundingbox.zmax+=gb.shiftz;
		// calculate distance
		double distanceb = communication_distance2_of_aabb_to_cell(boundingbox,node);
		if (distance > distanceb) distance = distanceb;
	}
	}
	}
	return distance;
}

extern double opening_angle2;
void communication_prepare_essential_cell_for_proc(struct cell* node, int proc){
	// Add essential cell to tree_essential_send
	if (tree_essential_send_N[proc]>=tree_essential_send_Nmax[proc]){
		tree_essential_send_Nmax[proc] += 32;
		tree_essential_send[proc] = realloc(tree_essential_send[proc],sizeof(struct cell)*tree_essential_send_Nmax[proc]);
	}
	tree_essential_send[proc][tree_essential_send_N[proc]] = (*node);
	tree_essential_send_N[proc]++;
		
	double distance2 = communication_distance2_of_proc_to_node(proc,node);
	if (((node->w*node->w)> opening_angle2*distance2)) {
		for (int o=0;o<8;o++){
			struct cell* d = node->oct[o];
			if (d==NULL) continue;
			communication_prepare_essential_cell_for_proc(d,proc);
		}
	}
}

void communication_prepare_essential_tree(struct cell* root){
	if (root==NULL) return;
	// Find out which cells are needed by every other node
	for (int i=0; i<mpi_num; i++){
		if (i==mpi_id) continue;
		communication_prepare_essential_cell_for_proc(root,i);	
	}
}

void communication_distribute_essential_tree(){
	// Distribute the number of cells to be transferred.
	for (int i=0;i<mpi_num;i++){
		MPI_Scatter(tree_essential_send_N, 1, MPI_INT, &(tree_essential_recv_N[i]), 1, MPI_INT, i, MPI_COMM_WORLD);
	}
	// Allocate memory for incoming tree_essential
	for (int i=0;i<mpi_num;i++){
		if  (i==mpi_id) continue;
		while (tree_essential_recv_Nmax[i]<tree_essential_recv_N[i]){
			tree_essential_recv_Nmax[i] += 32;
			tree_essential_recv[i] = realloc(tree_essential_recv[i],sizeof(struct cell)*tree_essential_recv_Nmax[i]);
		}
	}

	
	// Exchange tree_essential via MPI.
	// Using non-blocking receive call.
	MPI_Request request[mpi_num];
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (tree_essential_recv_N[i]==0) continue;
		MPI_Irecv(tree_essential_recv[i], tree_essential_recv_N[i], mpi_cell, i, i*mpi_num+mpi_id, MPI_COMM_WORLD, &(request[i]));
	}
	// Using blocking send call.
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (tree_essential_send_N[i]==0) continue;
		MPI_Send(tree_essential_send[i], tree_essential_send_N[i], mpi_cell, i, mpi_id*mpi_num+i, MPI_COMM_WORLD);
	}
	// Wait for all tree_essential to be received.
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (tree_essential_recv_N[i]==0) continue;
		MPI_Status status;
		MPI_Wait(&(request[i]), &status);
	}
	// Add tree_essential to local tree
	for (int i=0;i<mpi_num;i++){
		for (int j=0;j<tree_essential_recv_N[i];j++){
			tree_add_essential_node(&(tree_essential_recv[i][j]));
		}
	}
	// Bring everybody into sync, clean up. 
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<mpi_num;i++){
		tree_essential_send_N[i] = 0;
		tree_essential_recv_N[i] = 0;
	}
}
#endif


#endif // MPI
