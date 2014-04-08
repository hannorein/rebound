/**
 * @file 	collision_resolve.h
 * @brief 	Resolve a single collision. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	These functions resolve a single collision
 * of two colliding particles using momentum und energy conservation.
 * The coefficient resitution can be set with the variable 
 * coefficient_of_restitution or, when a velocity dependent 
 * coefficient of restitution is required with the function pointer
 * coefficient_of_restitution_for_velocity.
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
#ifndef _COLLISION_RESOLVE_H
#define _COLLISION_RESOLVE_H
#include "boundaries.h"
#include "collisions.h"

/**
 * Collision structure of one single collisions
 * Used to save a collision during collision search. 
 */
struct collision{
	int p1;			/**< First colliding particle. */
	int p2;			/**< Second colliding particle. */
	struct ghostbox gb;	/**< Ghostbox (of particle p1). */
#if defined(COLLISIONS_SWEEP) || defined(COLLISIONS_SWEEPPHI)
	double time;		/**< Time of collision. */
	int crossing;		/**< Collision occurs at the interface of two sweep boxes. */
#endif // COLLISIONS_SWEEP
	int ri;	 		/**< Index of rootcell (Needed for MPI). */
} collision;

extern double coefficient_of_restitution;	/**< Constant coefficient of restitution. Only used when coefficient_of_restitution_for_velocity not set by user. */
extern double minimum_collision_velocity;	/**< Minimal collision velocity. Needed to avoid particles slowly sinking into each other while constantly touching each other. */
extern double (*coefficient_of_restitution_for_velocity) (double); /**< Function pointer to a function that returns the coefficient of restitution as a function of velocity. */


/**
 * Resolve a single collision assuming a hardsphere collision model (no super-particle).
 * @param c Collision to resolve.
 */
void collision_resolve_hardsphere(struct collision c);	

/**
 * Function pointer to collision resolve function. 
 * This can be overwritten by the user to allow for implementation of super-particles
 * or anything where collisions are not resolved exactly.
 * @param c Collision to resolve.
 */
extern void (*collision_resolve) (struct collision);

#endif // _COLLISION_RESOLVE_H

