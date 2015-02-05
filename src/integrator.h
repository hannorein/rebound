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
#ifndef _INTEGRATOR_EULER_H
#define _INTEGRATOR_EULER_H
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
 * This parameter controls the accuracy of an adaptive integrator.
 * Default is 0 (non-adaptive).
 **/
extern double integrator_epsilon;

/*
 * The minimum timestep to be used in an adaptive integrator.
 * Default is 0 (no minimal timestep).
 **/
extern double integrator_min_dt;

#ifdef INTEGRATOR_IAS15 // MEGNO Routines are currently only implemented for IAS15
/* 
 * Init the MEGNO particles
 **/
void integrator_megno_init(double delta);

/*
 * Returns the current value of <Y>
 **/
double integrator_megno();

/*
 * Returns the largest Lyapunov characteristic number (LCN), or maximal Lyapunov exponent
 **/
double integrator_lyapunov();

#endif // INTEGRATOR_IAS15

#endif
