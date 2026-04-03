/**
 * @file    transformations.h
 * @brief   Transformations back and forth between various coordinate systems
 * @author  Hanno Rein <hanno@hanno-rein.de>
 *
 * @section LICENSE
 * Copyright (c) 2017 Hanno Rein, Dan Tamayo.
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

#ifndef _TRANSFORMATIONS_H
#define _TRANSFORMATIONS_H

// Jacobi
DLLEXPORT void reb_transformations_inertial_to_jacobi_posvel(const struct reb_particle* const particles, struct reb_particle* const p_j, const struct reb_particle* const p_mass, const size_t N, const size_t N_active); // p_mass: Should be the same particles array as ps for real particles. If passing variational particles in ps, p_mass should be the corresponding array of real particles.
DLLEXPORT void reb_transformations_inertial_to_jacobi_posvelacc(const struct reb_particle* const particles, struct reb_particle* const p_j, const struct reb_particle* const p_mass, const size_t N, const size_t N_active);
DLLEXPORT void reb_transformations_inertial_to_jacobi_acc(const struct reb_particle* const particles, struct reb_particle* const p_j,const struct reb_particle* const p_mass, const size_t N, const size_t N_active);
DLLEXPORT void reb_transformations_jacobi_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_j, const struct reb_particle* const p_mass, const size_t N, const size_t N_active);
DLLEXPORT void reb_transformations_jacobi_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_j, const struct reb_particle* const p_mass, const size_t N, const size_t N_active);
DLLEXPORT void reb_transformations_jacobi_to_inertial_acc(struct reb_particle* const particles, const struct reb_particle* const p_j, const struct reb_particle* const p_mass, const size_t N, const size_t N_active);


// Democratic heliocentric coordinates
DLLEXPORT void reb_transformations_inertial_to_democraticheliocentric_posvel(const struct reb_particle* const particles, struct reb_particle* const p_h, const size_t N, const size_t N_active);
DLLEXPORT void reb_transformations_democraticheliocentric_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_h, const size_t N, const size_t N_active);
DLLEXPORT void reb_transformations_democraticheliocentric_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_h, const size_t N, const size_t N_active);

// WHDS
DLLEXPORT void reb_transformations_inertial_to_whds_posvel(const struct reb_particle* const particles, struct reb_particle* const p_h, const size_t N, const size_t N_active);
DLLEXPORT void reb_transformations_whds_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_h, const size_t N, const size_t N_active);
DLLEXPORT void reb_transformations_whds_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_h, const size_t N, const size_t N_active);

// Barycentric coordinates
DLLEXPORT void reb_transformations_inertial_to_barycentric_posvel(const struct reb_particle* const particles, struct reb_particle* const p_b, const size_t N, const size_t N_active);
DLLEXPORT void reb_transformations_barycentric_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_b, const size_t N, const size_t N_active);
DLLEXPORT void reb_transformations_barycentric_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_b, const size_t N, const size_t N_active);
DLLEXPORT void reb_transformations_inertial_to_barycentric_acc(const struct reb_particle* const particles, struct reb_particle* const p_b, const size_t N, const size_t N_active);
DLLEXPORT void reb_transformations_barycentric_to_inertial_acc(struct reb_particle* const particles, const struct reb_particle* const p_b, const size_t N, const size_t N_active);
#endif
