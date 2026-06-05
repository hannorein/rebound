/**
 * integrator_ias15.h: The high order IAS15 integrator.
 * 
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
#ifndef _INTEGRATOR_IAS15_H
#define _INTEGRATOR_IAS15_H

extern const struct reb_integrator reb_integrator_ias15;

struct reb_integrator_ias15_state {
    double epsilon;                 // Precision control parameter
    double min_dt;                  // Minimal timestep
#define REB_INTEGRATOR_IAS15_ADAPTIVEMODE(X,Y) \
    X(Y, 0, INDIVIDUAL)     /* fractional error is calculated separately for each particle              */ \
    X(Y, 1, GLOBAL)         /* fractional error is calculated globally (was default until 01/2024)      */ \
    X(Y, 2, PRS23)          /* Pham, Rein & Spiegel (2023) timestep criterion (default since 01/2024)   */ \
    X(Y, 3, AARSETH85)       /* Aarseth (1985) timestep criterion                                        */
    enum {
        REB_GENERATE_ENUM(REB_INTEGRATOR_IAS15_ADAPTIVEMODE)
    } adaptive_mode;                    // Determines how the timestep is chosen
#define REB_INTEGRATOR_IAS15_OPTIMIZATION(X,Y) \
    X(Y, 0, DEFAULT)     /* Default */ \
    X(Y, 1, NO_CS)       /* No compensated summation */
    enum {
        REB_GENERATE_ENUM(REB_INTEGRATOR_IAS15_OPTIMIZATION)
    } optimization;                     // Optimization mode
    uint64_t iterations_max_exceeded;   // Counter how many times the iteration did not converge. 
    size_t N_allocated;          
    double* REB_RESTRICT at;
    double* REB_RESTRICT x0;
    double* REB_RESTRICT v0;
    double* REB_RESTRICT a0;
    double* REB_RESTRICT csx;
    double* REB_RESTRICT csv;
    double* REB_RESTRICT csa0;
    // The following are reb_dp7 pointers. See implementation for details.
    double* REB_RESTRICT g;
    double* REB_RESTRICT b;
    double* REB_RESTRICT csb;             // Compensated summation storage for b
    double* REB_RESTRICT e;
    double* REB_RESTRICT br;              // Used for resetting the b coefficients if a timestep gets rejected
    double* REB_RESTRICT er;              // Same for e coefficients
};

// Returns the gravitational timescale as calculated in Pham, Rein, Spiegel (2023). Useful for setting the initial IAS15 timestep.
REB_API double reb_integrator_ias15_timescale(struct reb_simulation* r);

#endif
