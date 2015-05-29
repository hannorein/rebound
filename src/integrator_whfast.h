/**
 * @file 	integrator_whfast.h
 * @brief 	Interface for numerical particle integrator
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Hanno Rein, Daniel Tamayo
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
#ifndef _INTEGRATOR_WHFAST_H
#define _INTEGRATOR_WHFAST_H
void integrator_whfast_part1(void);
void integrator_whfast_part2(void);
void integrator_whfast_synchronize(void);
void integrator_whfast_reset(void);

/*
 * This variable turns on/off various symplectic correctors.
 * 0 (default): turns off all correctors
 * 3: uses third order (two-stage) corrector 
 * 5: uses fifth order (four-stage) corrector 
 * 7: uses seventh order (six-stage) corrector 
 * 11: uses eleventh order (ten-stage) corrector 
 */
extern unsigned int integrator_whfast_corrector;

/* 
 * Setting this flag to one will recalculate Jacobi coordinates 
 * from the particle structure in the next timestep only. 
 * Then the flag gets set back to 0. If you want to change 
 * particles after every timestep, you also need to set this 
 * flag to 1 before every timestep.
 * Default is 0.
 **/ 
extern unsigned int integrator_whfast_particles_modified;

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
 * Note that whenever you change particles manually, you need to
 * set integrator_whfast_particles_modified to 1.
 **/
extern unsigned int integrator_whfast_synchronize_manually;
#endif
