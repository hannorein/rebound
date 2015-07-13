/**
 * @file 	collisions.h
 * @brief 	Collision search. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	The code supports different methods for the collision
 * detection. They all use this common interface. To detect and resolve 
 * collisions correctly, positions and velocities of the particles have 
 * to be synchronized in time. For that reason collisions_search() is
 * called at the end of the DKD timestep.  
 * 
 * 
 * @section LICENSE
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
#ifndef _COLLISIONS_H
#define _COLLISIONS_H
#include "boundary.h"
/**
 * Collision structure of one single collisions
 * Used to save a collision during collision search. 
 */
struct collision{
	int p1;			/**< First colliding particle. */
	int p2;			/**< Second colliding particle. */
	struct reb_ghostbox gb;	/**< Ghostbox (of particle p1). */
#if defined(COLLISIONS_SWEEP) || defined(COLLISIONS_SWEEPPHI)
	double time;		/**< Time of collision. */
	int crossing;		/**< Collision occurs at the interface of two sweep boxes. */
#endif // COLLISIONS_SWEEP
	int ri;	 		/**< Index of rootcell (Needed for MPI). */
} collision;

/**
 * Search for collisions and save them.
 * This routine searches for all collisions and stores them to 
 * resolve them when called collisions_resolve().
 */
void collisions_search(struct reb_context* const r);
/**
 * Resolve all collisions.
 * This function resolve all previously found collisions.
 */
void collisions_resolve(struct reb_context* const r);


/**
 * Just returns the constant coefficient of restitution in the REBOUND struct.
 */
double collisions_constant_coefficient_of_restitution_for_velocity(const struct reb_context* const r, double v);


/**
 * Resolve a single collision assuming a hardsphere collision model (no super-particle).
 * @param c Collision to resolve.
 */
void collision_resolve_hardsphere(struct reb_context* const r, struct collision c);	
#endif // _COLLISIONS_H
