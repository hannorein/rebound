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
/**
 * Search for collisions and save them.
 * This routine searches for all collisions and stores them to 
 * resolve them when called collisions_resolve().
 */
void collisions_search();
/**
 * Resolve all collisions.
 * This function resolve all previously found collisions.
 */
void collisions_resolve();

/**
 * Maximum radius of the particles in this simulation.
 * The value is set automatically in problem_init().
 */
extern double collisions_max_r; 

/**
 * Radius of the second largest particle in this simulation.
 * This radius is needed to decide if a tree cell needs to
 * be opened or not.
 * The value is set automatically in problem_init().
 */
extern double collisions_max2_r; 

/**
 * Particle between 0 and N_collisions (excluding) will not be searched for collisions and not shuffled around.
 */
extern int N_collisions;
#endif // _COLLISIONS_H
