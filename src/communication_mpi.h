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

void mpi_init(int argc, char** argv);
void add_particle_remote(struct particle pt, int proc_id);


#endif // MPI
#endif // _COMMUNICATION_MPI_H
