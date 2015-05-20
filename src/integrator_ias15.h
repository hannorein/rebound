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
void integrator_ias15_part1(void);
void integrator_ias15_part2(void);
void integrator_ias15_synchronize(void);
void integrator_ias15_reset(void);

/*
 * This parameter controls the accuracy of the integrator.
 * Set to 0 to make IAS15 a non-adaptive integrator.
 * Default: 1e-9.
 **/
extern double integrator_ias15_epsilon;

/*
 * The minimum timestep to be used in the adaptive integrator.
 * Default is 0 (no minimal timestep).
 **/
extern double integrator_ias15_min_dt;
#endif
