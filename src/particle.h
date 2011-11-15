/**
 * @file 	particle.h
 * @brief 	Particle structure and main particle routines.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
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
#ifndef _PARTICLE_H
#define _PARTICLE_H

#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
struct cell;
#endif // TREE

/**
 * Particle structure.
 * @details This structure is used to represent one particle. Additional particle
 * properties should be added here. Note that when the data structure is changed, 
 * one must also update the equivalent declaration for MPI in communications_mpi.c.
 */
struct particle {
	double x;	/**< x-position of the particle. */
	double y;	/**< y-position of the particle. */
	double z;	/**< z-position of the particle. */
	double vx;	/**< x-velocity of the particle. */
	double vy;	/**< y-velocity of the particle. */
	double vz;	/**< z-velocity of the particle. */
	double ax;	/**< x-acceleration of the particle. */
	double ay;	/**< y-acceleration of the particle. */
	double az;	/**< z-acceleration of the particle. */
	double m;	/**< Mass of the particle. */
#ifndef COLLISIONS_NONE
	double r; 	/**< Radius of the particle. */
	double lastcollision;	/**< Last time the particle had a physical collision. */
#endif // COLLISIONS_NONE
#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
	struct cell* c;		/**< Pointer to the cell the particle is currently in. */
#endif // TREE
};

/**
 * Main particle array.
 * This contains all particles on this node.
 */
extern struct particle* particles;

/** 
 * Adds a particle to the simulation. 
 * @details If a tree is used, it also adds the particle to the tree.
 *          If MPI is used, a particle that does not belong to the current node is 
 *          put into the corresponding send queue particles_send.
 * @param pt Particle to be added.
 */
void particles_add(struct particle pt);

/** 
 * Same as particles_add() but inserts particles at given position. 
 * @param pt Particle to be added.
 * @param pos New position.
 */
void particles_add_fixed(struct particle pt,int pos);

/**
 * Add a particle to the particle structure on the current node.
 * Do not distribute particles.
 * @param pt Particle to be added.
 */
void particles_add_local(struct particle pt);

/**
 * Returns the index of the rootbox for the current particles based on its position.
 * @return Index of the rootbox.
 * @param pt Particle to be checked.
 */
int particles_get_rootbox_for_particle(struct particle pt);

#endif // _PARTICLE_H
