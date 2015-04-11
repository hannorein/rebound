/**
 * @file 	integrator_mikkola.h
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
#ifndef _INTEGRATOR_MIKKOLA_H
#define _INTEGRATOR_MIKKOLA_H
void integrator_mikkola_part1();
void integrator_mikkola_part2();
void integrator_mikkola_synchronize();
void integrator_mikkola_reset();

/*
 * This variable turns on/off various symplectic correctors.
 * 0 (default): turns off all correctors
 * 3: uses third order (two-stage) corrector 
 * 5: uses fifth order (four-stage) corrector 
 * 7: uses seventh order (six-stage) corrector 
 */
extern unsigned int integrator_mikkola_corrector;

/* 
 * Flag determining if the integrator needs to recalculate the Jacobi
 * coordinates of each particle at every timestep. Set this to 1 if 
 * the masses of all particles stay constant during the entire simulation 
 * and the positions and velocities of particles are not changed by 
 * the user between timesteps.
 * Setting this to 1 results in a speed and accuracy increase.
 * Default is 0.
 **/ 
extern unsigned int integrator_mikkola_persistent_particles;

/* 
 * Flag overwriting the effect of integrator_mikkola_persistent_particle
 * for the next step only. Setting this flag to one will recalculate 
 * Jacobi coordinates from the particle structure in the next timestep
 * only. 
 * Default is 0.
 **/ 
extern unsigned int integrator_mikkola_particles_modified;

/*
 * Flag determining if the integrator produces synchronized outputs at
 * the end of the timestep. Setting this to 1 results in a speedup.
 * The general procedure to use this is:
 *   set integrator_synchronize_manually = 1
 *   run integrator for many steps.
 *   call integrator_synchronize()
 *   call output routines
 *   (continue with integration)
 *  
 * Default is 0 (produces synchronized outputs at every timestep).
 * Note that setting this flag to 1 implicitly also sets
 * integrator_mikkola_persistent_particles=1. 
 * This means you cannot change the particle mass/positions 
 * between timesteps unless you call integrator_synchronize();
 **/
extern unsigned int integrator_mikkola_synchronize_manually;
#endif
