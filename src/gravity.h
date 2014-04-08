/**
 * @file 	gravity.h
 * @brief 	Calculate gravitational forces. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	The code supports different methods for calculating the
 * gravitational forces. They all use this common interface. It is assumed
 * that the gravitational forces are independent of velocity. They are 
 * calculated by gravity_calculate_acceleration() which is called in the 
 * middle (K) part of the DKD timestepping scheme.  
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
#ifndef _GRAVITY_DIRECT_H
#define _GRAVITY_DIRECT_H
/**
  * The function loops over all ghostboxs and calls calculate_forces_for_particle() to sum up the forces on each particle.
  * Calculate all the gravitational acceleration for all particles.
  * Different methods implement this function in a different way.
  */
void gravity_calculate_acceleration();
#endif
