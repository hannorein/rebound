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
 * The first half of the integrator step.
 * This function is called at the beginning of the timestep. It 
 * advances the positions by 1/2 timestep.
 */
void integrator_part1();
/*
 * The second half of the integrator step.
 * This function is called after gravitational (and non-gravitational) 
 * forces for each particle have been calculated. It advances the 
 * velocity by 1 timestep and the positions by 1/2 timestep.
 * At the end of this function, the positions and velocities are in
 * sync which is needed for collision detection.
 */
void integrator_part2();
 

/* 
 * Flag determining if the integrator needs to consider velocity 
 * dependent forces. This is only relevant for IAS15.
 * Default is 1.
 **/ 
extern int integrator_force_is_velocitydependent;

/* 
 * Flag determining if the integrator needs to recalculate the Jacobi
 * coordinates of each particle at every timestep. This is only relevant 
 * for MIKKOLA as of now. Set this to 1 if the masses of all particles 
 * stay constant during the entire simulation and the positions and 
 * velocity of particles are not changed between timesteps.
 * Setting this to 1 results in a speed and accuracy increase.
 * Default is 0.
 **/ 
extern int integrator_intertial_frame;

/*
 * Flag determining if the integrator produces synchronized outputs at
 * the end of the timestep. Setting this to 0 results in a speedup.
 * The general procedure to use this is:
 *   set integrator_synchronized = 1
 *   run integrator_part1()
 *   set integrator_synchronized = 0
 *   run integrator_part2()
 *   run integrator_part1()
 *   (repeat last two steps many times until output is required)
 *   set integrator_synchronized = 1
 *   run integrator_part2()
 *   output
 *  
 * Default is 1 (produces synchronized outputs at every timestep).
 **/
extern unsigned int integrator_synchronized;

/*
 * This parameter controls the accuracy of an adaptive integrator.
 * Default is 0 (non-adaptive).
 **/
extern double integrator_epsilon;

/*
 * The minimum timestep to be used in an adaptive integrator.
 * Default is 0 (no minimal timestep).
 **/
extern double integrator_min_dt;

/* 
 * Cleanup all temporarily stored values.
 **/
void integrator_reset();

#endif
