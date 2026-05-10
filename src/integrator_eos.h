/**
 * integrator_eos.h: The Embedded Operator Splitting Integrator
 *
 * Copyright (c) 2019 Hanno Rein
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
#ifndef _INTEGRATOR_EOS_H
#define _INTEGRATOR_EOS_H

extern const struct reb_integrator reb_integrator_eos;

// Available methods for EOS Integrator
enum REB_INTEGRATOR_EOS_TYPE {
#define REB_INTEGRATOR_EOS_TYPE(X,Y) \
    X(Y, 0, LF)  \
    X(Y, 1, LF4) \
    X(Y, 2, LF6) \
    X(Y, 3, LF8) \
    X(Y, 4, LF4_2)   \
    X(Y, 5, LF8_6_4) \
    X(Y, 6, PLF7_6_4)\
    X(Y, 7, PMLF4)   \
    X(Y, 8, PMLF6)   
    REB_GENERATE_ENUM(REB_INTEGRATOR_EOS_TYPE)
};

struct reb_integrator_eos_state {
    enum REB_INTEGRATOR_EOS_TYPE phi0;         // Outer operator splitting method
    enum REB_INTEGRATOR_EOS_TYPE phi1;         // Inner operator splitting method
    unsigned int n;                 // Number of inner splittings per outer splitting
    unsigned int safe_mode;         // Combine Kick steps at beginning and end of timestep
};

#endif
