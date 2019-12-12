/**
 * @file 	gravity.h
 * @brief 	Calculate gravitational forces. 
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
#ifndef _GRAVITY_H
#define _GRAVITY_H
struct reb_simulation;

/**
  * The function loops over all ghostboxs and calls calculate_forces_for_particle() to sum up the forces on each particle.
  * Calculate all the gravitational acceleration for all particles.
  * Different methods implement this function in a different way.
  */
void reb_calculate_acceleration(struct reb_simulation* r);

/**
  * The function calculates the acceleration for the variational equations.
  */
void reb_calculate_acceleration_var(struct reb_simulation* r);


/**
  * The function calculates the jerk (derivative of the acceleration) and applies it to the particles' velocity.
  */
void reb_calculate_and_apply_jerk(struct reb_simulation* r, const double v);

#endif
