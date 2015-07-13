/**
 * @file 	boundary.h
 * @brief 	Handles different boundary conditions. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	The code supports different boundary conditions.
 * They all use this common interface. The function reb_boundary_check()
 * is called to check if particles have left the main box. If so, they are
 * shifted accordingly if the box is (shear) periodic and removed from the 
 * simulatiom.
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
#ifndef _BOUNDARIES_H
#define _BOUNDARIES_H

/**
 * This struct containes the relative position and velocity of a boundary box.
 * It is sometimes also used as the relative position and velocity of a 
 * particle to speed up calculation.
 */
struct reb_ghostbox{
	double shiftx;		/**< Relative x position */
	double shifty;		/**< Relative y position */
	double shiftz;		/**< Relative z position */
	double shiftvx;		/**< Relative x velocity */
	double shiftvy;		/**< Relative y velocity */
	double shiftvz;		/**< Relative z velocity */
};

/**
 * This function checks if any particle has left the main box.
 * If a particle left the box, it is shifted back in the box
 * for periodic boundary conditions or remove from the simulation.
 */
void reb_boundary_check(struct reb_context* r);

/**
 * Creates a ghostbox.
 * @param i Index in x direction.
 * @param j Index in y direction.
 * @param k Index in z direction.
 */
struct reb_ghostbox reb_boundary_get_ghostbox(struct reb_context* const r, int i, int j, int k);

/**
 * Return 1 if a particle is in the box, 0 otherwise.
 */
int reb_boundary_particle_is_in_box(const struct reb_context* const r, struct reb_particle p);

#endif
