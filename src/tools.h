/**
 * @file 	tools.h
 * @brief 	Tools for creating distributions.
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
#ifndef TOOLS_H
#define TOOLS_H
struct reb_simulation;



/**
 * @brief Returns deltad/delta 
 * @details Note, there is a typo in Gozdziewski et al 2001.
 * @param r REBOUND simulation to be considered.
 */
double reb_tools_megno_deltad_delta(struct reb_simulation* const r);

/**
 * @brief Update MEGNO after a successful timestep by adding dY (=ddelta/delta*dt)
 * @param r REBOUND simulation to be considered.
 * @param dY Increment for MEGNO Y
 */
void reb_tools_megno_update(struct reb_simulation* r, double dY);


/**
 * @brief Init random number generator based on time and process id.
 */
void reb_tools_init_srand(void);

#endif 	// TOOLS_H
