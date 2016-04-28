/**
 * @file 	integrator.h
 * @brief 	Interface for numerical particle integrators
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
#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H
struct reb_simulation;

/**
 * @brief The first half of the integrator step.
 * @details This function is called at the beginning of the timestep. It 
 * advances the positions by 1/2 timestep.
 */
void reb_integrator_part1(struct reb_simulation* r);

/**
 * @brief The second half of the integrator step.
 * @details This function is called after gravitational (and non-gravitational) 
 * forces for each particle have been calculated. It advances the 
 * velocity by 1 timestep and the positions by 1/2 timestep.
 * At the end of this function, the positions and velocities are in
 * sync which is needed for collision detection.
 */
void reb_integrator_part2(struct reb_simulation* r);


/** 
 * @brief This function updates the acceleration on all particles. 
 * @details It uses the current position and velocity data in the 
 * (struct reb_particle*) particles structure.
 */
void reb_update_acceleration(struct reb_simulation* r);

#endif
