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

#include "rebound.h"

DLLEXPORT void reb_integrator_whfast_from_inertial(struct reb_simulation* const r);   // Used by REBOUNDx
DLLEXPORT void reb_integrator_whfast_to_inertial(struct reb_simulation* const r); // Used by REBOUNDx
DLLEXPORT void reb_integrator_whfast_reset(struct reb_simulation* r);		// Used by REBOUNDx
DLLEXPORT int reb_integrator_whfast_init(struct reb_simulation* const r);    // Used by REBOUNDx, SABA
DLLEXPORT void reb_integrator_whfast_interaction_step(struct reb_simulation* const r, const double _dt);
DLLEXPORT void reb_integrator_whfast_jump_step(const struct reb_simulation* const r, const double _dt);
DLLEXPORT void reb_integrator_whfast_kepler_step(const struct reb_simulation* const r, const double _dt);
DLLEXPORT void reb_integrator_whfast_com_step(const struct reb_simulation* const r, const double _dt);
void reb_integrator_whfast_step(struct reb_simulation* r);
void reb_integrator_whfast_synchronize(struct reb_simulation* r);
void reb_integrator_whfast_calculate_jerk(struct reb_simulation* r);       ///< Calculates "jerk" term

#endif
