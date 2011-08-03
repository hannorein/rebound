/**
 * @file 	communication_mpi.h
 * @brief  	Handles communication between nodes using MPI.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	These routines handle the communication between
 * different nodes via the Message Passing Interface (MPI). 
 * There are two different types of communications
 * implemented at the moment:
 * - Distributing particles to the correct node.
 * - Creating, and distributing the essential tree to allow
 *   other nodes walk remote trees. Note that the opening 
 *   criteria is different for gravity and collision 
 *   tree walks.
 * 
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


#ifndef _COMMUNICATION_MPI_H
#define _COMMUNICATION_MPI_H
#ifdef MPI
#include "mpi.h"

int mpi_num;	/**< Total number of mpi nodes. */
int mpi_id;	/**< Unique ID of the current node in the range 0 <= mpid_id < mpi_num. */

/**
 * \defgroup mpistructures Data structures for MPI comminication.
 * These buffers are used for communicating with other nodes.
 * Each buffer is an array with the number of elements being equal
 * to the number of nodes. 
 * The number of particles/cells that have to be send to every other
 * node is saved in particles_send_N/tree_essential_send_N. The number
 * are different every time and have to be communicated to the other 
 * nodes before the actual communication of the particles/cells.
 * @{
 */
MPI_Datatype mpi_particle;				/**< MPI datatype corresponding to the C struct particle. */
extern struct particle** 	particles_send;		/**< Send buffer for particles. There is one buffer per node. */
extern int* 			particles_send_N;	/**< Current length of particle send buffer. */
extern int* 			particles_send_Nmax;	/**< Maximal length of particle send beffer before realloc() is needed. */
extern struct particle** 	particles_recv;		/**< Receive buffer for particles. There is one buffer per node. */
extern int* 			particles_recv_N;       /**< Current length of particle receive buffer. */
extern int* 			particles_recv_Nmax;    /**< Maximal length of particle receive beffer before realloc() is needed. */
#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
MPI_Datatype mpi_cell;						/**< MPI datatype corresponding to the C struct cell. */
struct cell;
extern struct cell** 		tree_essential_send;		/**< Send buffer for cells. There is one buffer per node. */
extern int* 			tree_essential_send_N;          /**< Current length of cell send buffer. */
extern int* 			tree_essential_send_Nmax;       /**< Maximal length of cell send beffer before realloc() is needed. */
extern struct cell** 		tree_essential_recv;            /**< Receive buffer for cells. There is one buffer per node. */
extern int* 			tree_essential_recv_N;          /**< Current length of cell receive buffer. */
extern int* 			tree_essential_recv_Nmax;       /**< Maximal length of cell receive beffer before realloc() is needed. */
#endif // TREE
/**@}*/

/**
 * Initializes MPI and sets up all necessary data structures.
 * @param argc Number of command line arguments. 
 * @param argc Command line arguments. 
 */
void communication_mpi_init(int argc, char** argv);

/**
 * Send particles in buffer particles_send to corresponding node. 
 * Receives particles from all nodes in buffer particles_recv and adds them
 * to the current simulation.
 */
void communication_mpi_distribute_particles();

/**
 * Places a particle in the send queue.  
 * @param pt Particle to be added to the send queue.
 * @param proc_id Particle will be send to this MPI node on next call of communication_mpi_distribute_particles();
 */
void communication_mpi_add_particle_to_send_queue(struct particle pt, int proc_id);

/**
 * Determine if the root box is local or if it is a copy of a remote node.
 * @param i Id of root box.
 */ 
int  communication_mpi_rootbox_is_local(int i);

#ifdef GRAVITY_TREE
/**
 * Send cells in buffer tree_essential_send to corresponding node. 
 * Receives cells from all nodes in buffer tree_essential_recv and adds them
 * to the non-local root boxes.
 */
void communication_mpi_distribute_essential_tree_for_gravity();

/**
 * Prepares the essential tree of a root box for communication with other nodes.
 * @param root The root cell under investigation.
 */
void communication_mpi_prepare_essential_tree_for_gravity(struct cell* root);
#endif // TREE

#ifdef COLLISIONS_TREE
/**
 * Send cells/particles in buffer tree_essential_send/particles_send to corresponding node. 
 * Receives cells/particles from all nodes in buffers. Does not insert particles 
 * into local tree.
 */
void communication_mpi_distribute_essential_tree_for_collisions();

/**
 * Prepares the essential tree/particles of a root box for communication with other nodes.
 * Adds copy of particles into particles_send.  
 * @param root The root cell under investigation.
 */
void communication_mpi_prepare_essential_tree_for_collisions(struct cell* root);
#endif //COLLISIONS_TREE


#endif // MPI
#endif // _COMMUNICATION_MPI_H
