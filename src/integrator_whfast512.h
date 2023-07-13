/**
 * @file 	integrator_whfast612.h
 * @brief 	Interface for numerical particle integrator
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *          Pejvak Javaheri <pejvak.javaheri@mail.utoronto.ca>
 * 
 * @section 	LICENSE
 * Copyright (c) 2023 Hanno Rein, Pejvak Javaheri
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
#ifndef _INTEGRATOR_WHFAST512_H
#define _INTEGRATOR_WHFAST512_H

#include "rebound.h"

void reb_integrator_whfast512_reset(struct reb_simulation* r);		
void reb_integrator_whfast512_part1(struct reb_simulation* r);		///< Internal function used to call a specific integrator
void reb_integrator_whfast512_part2(struct reb_simulation* r);		///< Internal function used to call a specific integrator
void reb_integrator_whfast512_synchronize(struct reb_simulation* r);	///< Internal function used to call a specific integrator
void reb_whfast512_kepler_solver(const struct reb_simulation* const r, struct reb_particle* const restrict p_j, const double M, unsigned int i, double _dt);   ///< Internal function (Main WHFast Kepler Solver)
void reb_whfast512_calculate_jerk(struct reb_simulation* r);       ///< Calculates "jerk" term

#endif
