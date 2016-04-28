/**
 * @file 	boundary.h
 * @brief 	Handles different boundary conditions. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
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
 * @brief This function checks if any particle has left the main box.
 * @details If a particle left the box, it is shifted back in the box
 * for periodic boundary conditions or remove from the simulation.
 * @param r REBOUND Simulation to consider
 */
void reb_boundary_check(struct reb_simulation* r);

/**
 * @brief Creates a ghostbox.
 * @param r REBOUND Simulation to consider
 * @param i Index in x direction.
 * @param j Index in y direction.
 * @param k Index in z direction.
 */
struct reb_ghostbox reb_boundary_get_ghostbox(struct reb_simulation* const r, int i, int j, int k);

/**
 * @details Return 1 if a particle is in the box, 0 otherwise.
 * @param r REBOUND Simulation to consider
 * @param p Particle to check
 */
int reb_boundary_particle_is_in_box(const struct reb_simulation* const r, struct reb_particle p);

#endif
