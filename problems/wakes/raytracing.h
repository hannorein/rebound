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
 * Calculates the opacity and the flux of reflexted light for a given geometry. 
 * @param N_rays Number of rays shot.
 * @param B elevation angle of observer (radians).
 * @param phi azimuthal angle of observer (radians).
 * @param sun_B elevation angle of light source (radians).
 * @param sun_phi azimuthal angle of light source (radians).
 * @param flux pointer to a double where the flux is written to
 * @param opacity pointer to a double where the opacity is written to
 */ 

void tree_raytrace(int N_rays, double B, double phi, double sun_B, double sun_phi, double* flux, double* opacity );


#endif // _RAYTRACING_H
