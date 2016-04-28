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


/**
 * Initializes MPI and sets up all necessary data structures.
 * @param argc Number of command line arguments. 
 * @param argc Command line arguments. 
 */
void reb_communication_mpi_init(struct reb_simulation* const r, int argc, char** argv);

/**
 * Send particles in buffer particles_send to corresponding node. 
 * Receives particles from all nodes in buffer particles_recv and adds them
 * to the current simulation.
 */
void reb_communication_mpi_distribute_particles(struct reb_simulation* const r);

/**
 * Places a particle in the send queue.  
 * @param pt reb_particle to be added to the send queue.
 * @param proc_id reb_particle will be send to this MPI node on next call of communication_mpi_distribute_particles();
 */
void reb_communication_mpi_add_particle_to_send_queue(struct reb_simulation* const r, struct reb_particle pt, int proc_id);

/**
 * Determine if the root box is local or if it is a copy of a remote node.
 * @param i Id of root box.
 */ 
int  reb_communication_mpi_rootbox_is_local(struct reb_simulation* const r, int i);

/**
 * Send cells in buffer tree_essential_send to corresponding node. 
 * Receives cells from all nodes in buffer tree_essential_recv and adds them
 * to the non-local root boxes.
 */
void reb_communication_mpi_distribute_essential_tree_for_gravity(struct reb_simulation* const r);

/**
 * Prepares the essential tree of a root box for communication with other nodes.
 * @param root The root cell under investigation.
 */
void reb_communication_mpi_prepare_essential_tree_for_gravity(struct reb_simulation* const r, struct reb_treecell* root);

/**
 * Send cells/particles in buffer tree_essential_send/particles_send to corresponding node. 
 * Receives cells/particles from all nodes in buffers. Does not insert particles 
 * into local tree.
 */
void reb_communication_mpi_distribute_essential_tree_for_collisions(struct reb_simulation* const r);

/**
 * Prepares the essential tree/particles of a root box for communication with other nodes.
 * Adds copy of particles into particles_send.  
 * @param root The root cell under investigation.
 */
void reb_communication_mpi_prepare_essential_tree_for_collisions(struct reb_simulation* const r, struct reb_treecell* root);


#endif // MPI
#endif // _COMMUNICATION_MPI_H
