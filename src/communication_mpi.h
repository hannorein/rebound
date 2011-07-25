#ifndef _COMMUNICATION_MPI_H
#define _COMMUNICATION_MPI_H
#ifdef MPI
#include "mpi.h"

int mpi_num;
int mpi_id;
MPI_Datatype mpi_particle;
#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
MPI_Datatype mpi_cell;
#endif

void communication_mpi_init(int argc, char** argv);
void communication_mpi_distribute_particles();
void communication_mpi_add_particle_to_send_queue(struct particle pt, int proc_id);
int  communication_mpi_rootbox_is_local(int i);

#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
struct cell;
void communication_mpi_distribute_essential_tree();
void communication_mpi_prepare_essential_tree(struct cell* root);
#endif 


#endif // MPI
#endif // _COMMUNICATION_MPI_H
