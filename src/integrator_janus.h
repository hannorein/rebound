/**
 * integrator_janus.h: Bit-wise reversible symplectic integrator
 * 
 * Copyright (c) 2017 Hanno Rein, Daniel Tamayo
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
#ifndef _INTEGRATOR_JANUS_H
#define _INTEGRATOR_JANUS_H

extern const struct reb_integrator reb_integrator_janus;

#define REB_PARTICLE_INT_TYPE int64_t
struct reb_particle_int {
    REB_PARTICLE_INT_TYPE x;
    REB_PARTICLE_INT_TYPE y;
    REB_PARTICLE_INT_TYPE z;
    REB_PARTICLE_INT_TYPE vx;
    REB_PARTICLE_INT_TYPE vy;
    REB_PARTICLE_INT_TYPE vz;
};

struct reb_integrator_janus_state {
    double scale_pos;       // Scale of position grid. Default 1e-16
    double scale_vel;       // Scale of velocity grid. Default 1e-16
    unsigned int order;     // Order: 2 (default), 4, 6, 8, 10 
    unsigned int recalculate_integer_coordinates_this_timestep;  // Set to 1 if particles have been modified

    // Internal use
    struct reb_particle_int* REB_RESTRICT p_int;
    size_t N_allocated;
};


#endif
