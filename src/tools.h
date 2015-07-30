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
#include "particle.h"
struct reb_simulation;



/**
 * @brief Initialize a particle on an orbit in the xy plane.
 * @param G Gravitational constant.
 * @param M Mass of the central object.
 * @param m Mass of the particle.
 * @param a Semi-major axis of the particle.
 * @param e Eccentricity of the particle.
 * @param omega Pericenter of the particle.
 * @param f true anomaly of the particle.
 * @return Returns a particle structure with the given orbital parameters. All other particle properties are 0 by default.
 */
struct reb_particle reb_tools_init_orbit2d(double G, double M, double m, double a, double e, double omega, double f);

/**
 * @brief Initialize a particle on a 3D orbit.  See Fig. 2.13 of Murray & Dermott Solar System Dynamics for diagram.
 * @param G Gravitational constant.
 * @param M Mass of the central object.
 * @param m Mass of the particle.
 * @param a Semi-major axis of the particle.
 * @param e Eccentricity of the particle.
 * @param i inclination of the particle to the reference plane.
 * @param Omega Longitude of the ascending node of the particle.
 * @param omega argument of pericenter of the particle.
 * @param f true anomaly of the particle.
 * @return Returns a particle structure with the given orbital parameters. All other particle properties are 0 by default.
 */
struct reb_particle reb_tools_init_orbit3d(double G, double M, double m, double a, double e, double i, double Omega, double omega, double f);

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
void reb_tools_init_srand();

#endif 	// TOOLS_H
