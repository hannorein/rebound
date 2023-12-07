/**
 * @file 	integrator_trace.h
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
#ifndef _INTEGRATOR_TRACE_H
#define _INTEGRATOR_TRACE_H
void reb_integrator_trace_part1(struct reb_simulation* r);          ///< Internal function used to call a specific integrator
void reb_integrator_trace_part2(struct reb_simulation* r);          ///< Internal function used to call a specific integrator
void reb_integrator_trace_synchronize(struct reb_simulation* r);    ///< Internal function used to call a specific integrator
void reb_integrator_trace_reset(struct reb_simulation* r);          ///< Internal function used to call a specific integrator
void reb_integrator_trace_inertial_to_dh(struct reb_simulation* r); ///< Internal in-place coordinate transformation
void reb_integrator_trace_dh_to_inertial(struct reb_simulation* r); ///< Internal in-place coordinate transformation

// Switching functions
double reb_integrator_trace_peri_switch_default(struct reb_simulation* const r, const unsigned int j);
double reb_integrator_trace_switch_fdot_peri(struct reb_simulation* const r, const unsigned int j);
double reb_integrator_trace_switch_vdiff_peri(struct reb_simulation* const r, const unsigned int j);

double reb_integrator_trace_switch_velocity(struct reb_simulation* const r, const unsigned int i, const unsigned int j);
#endif
