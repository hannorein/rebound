/**
 * integrator_whfast.h: Symplectic Wisdom-Holman integrator WHFast
 * 
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

extern const struct reb_integrator reb_integrator_whfast;

// WHFast Integrator (Rein & Tamayo 2015)
struct reb_integrator_whfast_state {
    unsigned int corrector;                                     // Order of first symplectic corrector: 0 (default - no corrector), 3, 5, 7, 11, 17.  
    unsigned int corrector2;                                    // 0: no second corrector, 1: use second corrector
#define REB_INTEGRATOR_WHFAST_KERNEL(X,Y) \
    X(Y, 0, DEFAULT) \
    X(Y, 1, MODIFIEDKICK) \
    X(Y, 2, COMPOSITION) \
    X(Y, 3, LAZY)
    enum {
        REB_GENERATE_ENUM(REB_INTEGRATOR_WHFAST_KERNEL)
    } kernel;                                                   // Kernel type. See Rein, Tamayo & Brown 2019 for details.                            
#define REB_INTEGRATOR_WHFAST_COORDINATES(X,Y) \
    X(Y, 0, JACOBI)                                         /* Jacobi coordinates (default)                   */ \
    X(Y, 1, DEMOCRATICHELIOCENTRIC)                         /* Democratic Heliocentric coordinates            */ \
    X(Y, 2, WHDS)                                           /* WHDS coordinates (Hernandez and Dehnen, 2017)  */ \
    X(Y, 3, BARYCENTRIC)                                    /* Barycentric coordinates                        */ 
    enum {
        REB_GENERATE_ENUM(REB_INTEGRATOR_WHFAST_COORDINATES)
    } coordinates;                                              // Coordinate system used in Hamiltonian splitting
    unsigned int safe_mode;                                     // 0: Drift Kick Drift scheme (default), 1: combine first and last sub-step.
    unsigned int keep_unsynchronized;                           // 1: continue from unsynchronized state after synchronization 

    // Internal use
    size_t N_allocated;
    struct reb_particle* REB_RESTRICT p_jh;     // Jacobi/heliocentric/WHDS coordinates
    size_t N_allocated_var;
    struct reb_particle* REB_RESTRICT p_jh_var; // Jacobi coordinates for variational equations
    size_t N_allocated_temp;
    struct reb_particle* REB_RESTRICT p_temp;   // Used for lazy implementer's kernel 
    unsigned int recalculate_coordinates_but_not_synchronized_warning;
};


REB_API void reb_integrator_whfast_from_inertial(struct reb_simulation* const r, struct reb_particle* p_jh, int coordinates);
REB_API void reb_integrator_whfast_to_inertial(struct reb_simulation* const r, struct reb_particle* p_jh, int coordinates);
REB_API int reb_integrator_whfast_init(struct reb_simulation* const r, struct reb_integrator_whfast_state* whfast);    // Used by REBOUNDx
REB_API void reb_integrator_whfast_interaction_step(struct reb_simulation* const r, struct reb_particle* p_jh, int coordinates, const double _dt);
REB_API void reb_integrator_whfast_jump_step(const struct reb_simulation* const r, struct reb_integrator_whfast_state* whfast, const double _dt); // Used by REBOUNDx
REB_API void reb_integrator_whfast_kepler_step(const struct reb_simulation* const r, struct reb_particle* p_jh, int coordinates, const double _dt);
REB_API void reb_integrator_whfast_com_step(const struct reb_simulation* const r, struct reb_particle* p_jh, const double _dt);
void reb_integrator_whfast_calculate_jerk(struct reb_simulation* r, struct reb_particle* jerk);

#endif
