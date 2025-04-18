/**
 * @file 	integrator_brace.h
 * @brief 	Interface for numerical particle integrator
 * @author 	Hanno Rein
 *
 * @section 	LICENSE
 * Copyright (c) 2023 Tiger Lu
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
#ifndef _INTEGRATOR_BRACE_H
#define _INTEGRATOR_BRACE_H
void reb_integrator_brace_part1(struct reb_simulation* r);          ///< Internal function used to call a specific integrator
void reb_integrator_brace_part2(struct reb_simulation* r);          ///< Internal function used to call a specific integrator
void reb_integrator_brace_synchronize(struct reb_simulation* r);    ///< Internal function used to call a specific integrator
void reb_integrator_brace_reset(struct reb_simulation* r);          ///< Internal function used to call a specific integrator
void reb_integrator_brace_dh_to_inertial(struct reb_simulation* r); ///< Ionternal in-place coordinate transformation
void reb_integrator_brace_inertial_to_dh(struct reb_simulation* r); ///< Internal in-place coordinate transformation
#endif
