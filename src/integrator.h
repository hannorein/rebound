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
 * @brief This function is used to initialize constants in some integrators. 
 * @details The function doesn't need to be called. Integrators will call it
 * from within the normal reb_integrator_part1() function. It is sometimes
 * called before an integration step is performed to ensure variables are 
 * set before a binary file is outputted.
 */
void reb_integrator_init(struct reb_simulation* r);
#endif
