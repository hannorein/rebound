/**
 * integrator_whfast512.h: The AVX512 accelerated symplectic integrator WHFast512 in ASM
 * 
 * Copyright (c) 2026 Hanno Rein, Rishit Dagli, Pejvak Javaheri
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
#ifndef _INTEGRATOR_WHFAST512_H
#define _INTEGRATOR_WHFAST512_H

#include "rebound.h"

extern const struct reb_integrator reb_integrator_whfast512;

struct reb_integrator_whfast512_state {
    unsigned int gr_potential;
    unsigned int N_systems;
    unsigned int corrector;             
    unsigned int concatenate_steps;

    // Internal use
    size_t N_allocated;
    double last_synchronization;
    void* data;
};

#endif
