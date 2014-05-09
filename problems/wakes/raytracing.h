/**
 * @file 	raytracing.c
 * @brief 	Tree-based ray tracing algorithm.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2014 Hanno Rein
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
#ifndef _RAYTRACING_H
#define _RAYTRACING_H

/**
 * Calculates the flux f. 
 * @return Returns 1 if simulation was restarted. 0 otherwise.
 * @param argc Number of command line arguments.
 * @param argv Array of command line arguments.
 */ 

void tree_get_transparency(double B, double phi, double sun_B, double sun_phi);


#endif _RAYTRACING_H
