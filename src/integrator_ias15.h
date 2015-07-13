/**
 * @file 	integrator_ias15.h
 * @brief 	Interface for numerical particle integrator
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Hanno Rein
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
#ifndef _INTEGRATOR_IAS15_H
#define _INTEGRATOR_IAS15_H
struct reb_context;
struct reb_particle; 

void reb_integrator_ias15_part1(struct reb_context* r);
void reb_integrator_ias15_part2(struct reb_context* r);
void reb_integrator_ias15_synchronize(struct reb_context* r);
void reb_integrator_ias15_reset(struct reb_context* r);

struct reb_context_integrator_ias15 {
	/**
	 * This parameter controls the accuracy of the integrator.
	 * Set to 0 to make IAS15 a non-adaptive integrator.
	 * Default: 1e-9.
	 **/
	double epsilon;

	/**
	 * The minimum timestep to be used in the adaptive integrator.
	 * Default is 0 (no minimal timestep).
	 **/
	double min_dt;
	
	/** 
	 * If 1: estimate the fractional error by max(acceleration_error)/max(acceleration), where max is take over all particles.
	 * If 0: estimate the fractional error by max(acceleration_error/acceleration).
	 **/
	unsigned int epsilon_global;


	// Internal data structures below. Nothing to be changed by the user.
	
	/**
	 * Count how many times the iteration did not converge. 
	 **/
	unsigned long iterations_max_exceeded;



	int N3allocated; 			// Size of allocated arrays.

	double* at;				// Temporary buffer for acceleration
	double* x0;				// Temporary buffer for position (used for initial values at h=0) 
	double* v0;				//                      velocity
	double* a0;				//                      acceleration
	double* csx;				//                      compensated summation
	double* csv;				//                      compensated summation

	double* g[7];
	double* b[7];
	double* e[7];

	// The following values are used for resetting the b and e coefficients if a timestep gets rejected
	double* br[7];
	double* er[7];
	double dt_last_success;			// Last accepted timestep (corresponding to br and er)

};
#endif
