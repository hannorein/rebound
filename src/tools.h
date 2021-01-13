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

#include <stdint.h>

struct reb_simulation;
struct reb_particles;

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
void reb_tools_init_srand(struct reb_simulation* r);

/**
 * @brief Convert angles for orbit routines
 */
double reb_tools_M_to_E(double e, double M);

/**
 * @brief Convert angles for orbit routines
 */
double reb_tools_M_to_f(double e, double M);

/**
 * @brief Kepler solver in Pal coordinates
 */
void reb_tools_solve_kepler_pal(double h, double k, double lambda, double* p, double* q);

/**
 * @brief Convert particle to Pal coordinates
 */
void reb_tools_particle_to_pal(double G, struct reb_particle p, struct reb_particle primary, double *a, double* lambda, double* k, double* h, double* ix, double* iy);

/**
 * @brief internal function to handle outputs for the Fast Simulation Restarter.
 */
void reb_fsr_heartbeat(struct reb_simulation* const r);

/**
 * @brief This function multiplies positions and velocities with a scalar. Note: modifies r.
 */
void reb_simulation_imul(struct reb_simulation* r, double scalar_pos, double scalar_vel);

/**
 * @brief This function adds the positions and velocities of particles in simulation r2 to those in simulation r. Note: modifies r.
 */
int reb_simulation_iadd(struct reb_simulation* r, struct reb_simulation* r2);

/**
 * @brief This function subtracts the positions and velocities of particles in simulation r2 from those in simulation r. Note: modifies r.
 */
int reb_simulation_isub(struct reb_simulation* r, struct reb_simulation* r2);

#endif 	// TOOLS_H
