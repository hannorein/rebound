/**
 * integrator_saba.h: SABA integrator family (Laskar and Robutel 2001)
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
#ifndef _INTEGRATOR_SABA_H
#define _INTEGRATOR_SABA_H

extern const struct reb_integrator reb_integrator_saba;

struct reb_integrator_saba_state {
#define REB_INTEGRATOR_SABA_TYPE(X,Y) \
    X(Y, 0, DEFAULT) \
    X(Y, 0x0,    1)         /* WH                                 */ \
    X(Y, 0x1,    2)         /* SABA2                              */ \
    X(Y, 0x2,    3)         /* SABA3                              */ \
    X(Y, 0x3,    4)         /* SABA4                              */ \
    X(Y, 0x100,  CM_1)      /* SABACM1 (Modified kick corrector)  */ \
    X(Y, 0x101,  CM_2)      /* SABACM2 (Modified kick corrector)  */ \
    X(Y, 0x102,  CM_3)      /* SABACM3 (Modified kick corrector)  */ \
    X(Y, 0x103,  CM_4)      /* SABACM4 (Modified kick corrector)  */ \
    X(Y, 0x200,  CL_1)      /* SABACL1 (lazy corrector)           */ \
    X(Y, 0x201,  CL_2)      /* SABACL2 (lazy corrector)           */ \
    X(Y, 0x202,  CL_3)      /* SABACL3 (lazy corrector)           */ \
    X(Y, 0x203,  CL_4)      /* SABACL4 (lazy corrector)           */ \
    X(Y, 0x4,    10_4)      /* SABA(10,4), 7 stages               */ \
    X(Y, 0x5,    8_6_4)     /* SABA(8,6,4), 7 stages              */ \
    X(Y, 0x6,    10_6_4)    /*SABA(10,6,4), 8 stages, default     */ \
    X(Y, 0x7,    H_8_4_4)   /*SABAH(8,4,4), 6 stages              */ \
    X(Y, 0x8,    H_8_6_4)   /*SABAH(8,6,4), 8 stages              */ \
    X(Y, 0x9,    H_10_6_4)  /*SABAH(10,6,4), 9 stages             */
    enum {
        REB_GENERATE_ENUM(REB_INTEGRATOR_SABA_TYPE)
    } type;                             // Type of integrator
    unsigned int safe_mode;             // Combine first and last sub-step
    unsigned int keep_unsynchronized;   // 1: continue from unsynchronized state after synchronization

    // Internal use
    size_t N_allocated;
    struct reb_particle* REB_RESTRICT p_jh;     // Jacobi/heliocentric/WHDS coordinates
    size_t N_allocated_temp;
    struct reb_particle* REB_RESTRICT p_temp;
};


#endif
