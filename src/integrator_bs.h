/**
 * integrator_bs.h: Bulirsch Stoer Integrator
 * 
 * Copyright (c) 2021 Hanno Rein
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
#ifndef _INTEGRATOR_BS_H
#define _INTEGRATOR_BS_H

extern const struct reb_integrator reb_integrator_bs;

struct reb_integrator_bs_state {
    double eps_abs; // Allowed absolute scalar error.
    double eps_rel; // Allowed relative scalar error.
    double min_dt;  // Minimum timestep
    double max_dt;  // Maximum teimstep

    // Internal use
    struct reb_ode* nbody_ode;  // ODE corresponding to N-body system
    int* sequence;              // stepsize sequence
    int* cost_per_step;         // overall cost of applying step reduction up to iteration k + 1, in number of calls.
    double* cost_per_time_unit; // cost per unit step.
    double* optimal_step;       // optimal steps for each order. 
    double* coeff;              // extrapolation coefficients.
    double dt_proposed;
    int first_or_last_step;
    int previous_rejected;
    int target_iter;
    int user_ode_needs_nbody;   // Do not set manually. Use needs_nbody in reb_ode instead.
};

int reb_integrator_bs_step_odes(struct reb_simulation* r, struct reb_integrator_bs_state* bs, double dt);
void reb_integrator_bs_update_particles(struct reb_simulation* r, const double* y);
void reb_integrator_bs_nbody_derivatives(struct reb_ode* ode, double* const yDot, const double* const y, double const t);
#endif
