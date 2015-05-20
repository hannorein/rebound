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

/*
 * Available integrator.
 */
typedef enum {
	IAS15 = 0,
	WHFAST = 1,
	SEI = 2,
	WH = 3,
	LEAPFROG = 4,
	HYBRID = 5,
	NONE = 6,
	} integrator_t;
/*
 * Variable setting the current integrator.
 */
extern integrator_t integrator;


/*
 * The first half of the integrator step.
 * This function is called at the beginning of the timestep. It 
 * advances the positions by 1/2 timestep.
 */
void integrator_part1(void);
/*
 * The second half of the integrator step.
 * This function is called after gravitational (and non-gravitational) 
 * forces for each particle have been calculated. It advances the 
 * velocity by 1 timestep and the positions by 1/2 timestep.
 * At the end of this function, the positions and velocities are in
 * sync which is needed for collision detection.
 */
void integrator_part2(void);
 

/* 
 * Flag determining if the integrator needs to consider velocity 
 * dependent forces. 
 * Default is 0.
 **/ 
extern unsigned int integrator_force_is_velocitydependent;


/*
 * Synchronize particles manually at end of timestep.
 */
void integrator_synchronize(void);

/* 
 * Cleanup all temporarily stored values.
 **/
void integrator_reset(void);

/* This function updates the acceleration on all particles. 
 * It uses the current position and velocity data in the 
 * (struct particle*) particles structure.
 * Note: this does currently not work with MPI or any TREE module.
 */
void integrator_update_acceleration(void);

#endif
