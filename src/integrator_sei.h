/**
 * @file 	integrator_sei.h
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
#ifndef _INTEGRATOR_SEI_H
#define _INTEGRATOR_SEI_H
void reb_integrator_sei_part1(struct reb_simulation* r);
void reb_integrator_sei_part2(struct reb_simulation* r);
void reb_integrator_sei_synchronize(struct reb_simulation* r);
void reb_integrator_sei_reset(struct reb_simulation* r);

struct reb_simulation_integrator_sei {
	double OMEGA;		/**< Epicyclic/orbital frequency.  */
	double OMEGAZ; 		/**< Epicyclic frequency in vertical direction. */

	double lastdt;		/**< Cached sin(), tan() for this value of dt.*/
	// Cache sin() tan() values.
	double sindt;
	double tandt;
	double sindtz;
	double tandtz;
};
#endif
