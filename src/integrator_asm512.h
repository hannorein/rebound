/**
 * integrator_asm512.h: The AVX512 accelerated symplectic integrator WHFast512 in ASM
 * 
 * Copyright (c) 2025 Hanno Rein
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
#ifndef _INTEGRATOR_ASM512_H
#define _INTEGRATOR_ASM512_H

#include "rebound.h"

extern const struct reb_integrator reb_integrator_asm512;

struct reb_integrator_asm512_state {
    unsigned int gr_potential;          // 1: Turn on GR potential of central object, 0 (default): no GR potential
    unsigned int N_systems;             // Number of systems to be integrator in parallel: 1 (default, up to 8 planets), 2 (up to 4 planets each), 4 (2 planets each)
    unsigned int corrector;
    unsigned int concatenate_steps;

    // Internal use
    size_t N_allocated;
    double last_synchronization;        // Time of last synchronization (required to advance com)
    void* data; // alligned SIMD data
};

#endif
