/**
 * integrator_mercurius.h: They hybrid symplectic Mercurius Integrator
 * 
 * Copyright (c) 2017 Hanno Rein
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
#ifndef _INTEGRATOR_MERCURIUS_H
#define _INTEGRATOR_MERCURIUS_H

extern const struct reb_integrator reb_integrator_mercurius;

struct reb_integrator_mercurius_state {
    double (*L) (const struct reb_simulation* const r, double d, double dcrit); // Switching function (default same as Mercury) 
    double r_crit_hill;                                 // Critical switching distance in units of Hill radii
    unsigned int safe_mode;                             // Combine Kick steps at beginning and end of timestep

    // Internal use
    enum {
        REB_INTEGRATOR_MERCURIUS_MODE_WH = 0,
        REB_INTEGRATOR_MERCURIUS_MODE_ENCOUNTER = 1,
    } mode;
    size_t encounter_N;             // Number of particles currently having an encounter
    size_t encounter_N_active;      // Number of active particles currently having an encounter
    unsigned int tponly_encounter;  // 0 if any encounters are between two massive bodies. 1 if encounters only involve test particles
    size_t N_allocated;
    size_t N_allocated_additional_forces;
    size_t N_allocated_dcrit;       // Current size of dcrit arrays
    double* dcrit;                  // Precalculated switching radii for particles
    struct reb_particle* REB_RESTRICT particles_backup; //  contains coordinates before Kepler step for encounter prediction
    struct reb_particle* REB_RESTRICT particles_backup_additional_forces; // contains coordinates before Kepler step for encounter prediction
    size_t* encounter_map;          // Map to represent which particles are integrated with ias15
    struct reb_vec3d com_pos;       // Used to keep track of the center of mass during the timestep
    struct reb_vec3d com_vel;
};

#endif
