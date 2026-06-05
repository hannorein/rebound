/**
 * integrator_trace.h: TRACE integrator
 *
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

extern const struct reb_integrator reb_integrator_trace;

struct reb_integrator_trace_state {
    int (*S) (struct reb_simulation* const r, const size_t i, const size_t j);
    int (*S_peri) (struct reb_simulation* const r, const size_t j);

#define REB_INTEGRATOR_TRACE_PERIMODE(X,Y) \
    X(Y, 0, PARTIAL_BS) \
    X(Y, 1, FULL_BS) \
    X(Y, 2, FULL_IAS15)
    enum {
        REB_GENERATE_ENUM(REB_INTEGRATOR_TRACE_PERIMODE)
    } peri_mode;

#define REB_INTEGRATOR_TRACE_COORDINATES(X,Y) \
    X(Y, 1, DEMOCRATICHELIOCENTRIC) \
    X(Y, 2, WIDEBINARY)
    enum {
        REB_GENERATE_ENUM(REB_INTEGRATOR_TRACE_COORDINATES)
    } coordinates;

    double r_crit_hill;
    double peri_crit_eta;
    double r_crit_WB; // critical distance (units of hill radius) for close encounters with the binary in WB coordinates

    // Internal use
    enum {
        REB_INTEGRATOR_TRACE_MODE_INTERACTION = 0, // Interaction step
        REB_INTEGRATOR_TRACE_MODE_KEPLER = 1,      // Kepler step
        REB_INTEGRATOR_TRACE_MODE_FULL = 3,        // Doing everything in one step
    } mode;
    size_t encounter_N;                 // Number of particles currently having an encounter
    size_t encounter_N_active;          // Number of active particles currently having an encounter

    size_t N_allocated;
    size_t N_allocated_additional_forces;
    unsigned int tponly_encounter;      // 0 if any encounters are between two massive bodies. 1 if encounters only involve test particles

    struct reb_particle* REB_RESTRICT particles_backup; //  Contains coordinates before the entire step
    struct reb_particle* REB_RESTRICT particles_backup_kepler; //  Contains coordinates before kepler step
    struct reb_particle* REB_RESTRICT particles_backup_additional_forces; // For additional forces

    size_t* encounter_map;              // Map to represent which particles are integrated with BS
    size_t* encounter_map_backup;       // Contains encounter map from after pre-ts check. Used to retain memory of CEs flagged at this step.
    struct reb_vec3d com_pos;           // Used to keep track of the centre of mass during the timestep
    struct reb_vec3d com_vel;

    int* current_Ks; // Tracking K_ij for the entire timestep
    unsigned int current_C; // Tracking C for the entire timestep
    unsigned int force_accept; // Force accept for irreversible steps: collisions and adding particles
};

// Built in trace switching functions

REB_API int reb_integrator_trace_switch_peri_default(struct reb_simulation* const r, const size_t j);
REB_API int reb_integrator_trace_switch_peri_none(struct reb_simulation* const r, const size_t j);
REB_API int reb_integrator_trace_switch_default(struct reb_simulation* const r, const size_t i, const size_t j);

#endif
