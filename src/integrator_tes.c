/**
 * @file 	integrator_tes.c
 * @brief 	Terrestrial Exoplanet Simulator (TES) integration scheme
 * @author 	Peter Bartram <p.bartram@soton.ac.uk>
 * 
 * @section 	LICENSE
 * Copyright (c) 2021 Peter Bartram
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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "rebound.h"
#include "integrator_tes.h"
#include "dhem.h"
#include "radau.h"
#include "radau_step.h"
#include "UniversalVars.h"

void reb_integrator_tes_part1(struct reb_simulation* r){} // unused

void reb_integrator_tes_part2(struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    uint32_t N = r->N;

    if(r->ri_tes.allocated_N != N)
    {
        r->ri_tes.allocated_N = N;
        struct reb_particle* const particles = r->particles;
        
        // Adding new mallocs for data that isnt under the sim data structure here.
        r->ri_tes.particles_dh = (struct particles*)malloc(sizeof(struct reb_particle)*r->N);


        r->ri_tes.sim = Simulation_Init(r, N);

        // Convert from inertial to dh coords.
        reb_move_to_com(r); // This can be removed once all inertial frames are enabled.
        reb_transformations_inertial_to_democraticheliocentric_posvel(particles, r->ri_tes.particles_dh, r->N, r->N);

        for(uint32_t i=1;i<N;i++) 
        {
            r->ri_tes.sim->mass[i] =     r->ri_tes.particles_dh[i].m;
            r->ri_tes.sim->Q_dh[3*i] =   r->ri_tes.particles_dh[i].x;
            r->ri_tes.sim->Q_dh[3*i+1] = r->ri_tes.particles_dh[i].y;
            r->ri_tes.sim->Q_dh[3*i+2] = r->ri_tes.particles_dh[i].z;
            r->ri_tes.sim->P_dh[3*i] =   r->ri_tes.particles_dh[i].vx*r->ri_tes.particles_dh[i].m;
            r->ri_tes.sim->P_dh[3*i+1] = r->ri_tes.particles_dh[i].vy*r->ri_tes.particles_dh[i].m;
            r->ri_tes.sim->P_dh[3*i+2] = r->ri_tes.particles_dh[i].vz*r->ri_tes.particles_dh[i].m; 
        }
        // Need this until above loop is changed to start at zero.
        r->ri_tes.sim->mass[0] = particles[0].m;

        // Default value - this need changing to be proper
        double t0 = 0;

        // Store configuration values in the simulation. (should be able to remove this section entirely eventually)
        r->ri_tes.sim->t0 = t0;

        UniversalVars_Init(r);
        dhem_Init(r, r->ri_tes.sim, r->ri_tes.orbital_period/r->ri_tes.recti_per_orbit, 9);
        dhem_InitialiseOsculatingOrbits(r, r->ri_tes.sim->Q_dh, r->ri_tes.sim->P_dh, r->ri_tes.sim->t0);
        Radau_Init(r);  
    }
    double dt_new = Radau_SingleStep(r, r->t, r->dt, r->dt_last_done);

    // update timestep
	r->t+=r->dt;
	r->dt_last_done = r->dt;
    r->dt = dt_new;
}

void reb_integrator_tes_synchronize(struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    uint32_t N = r->N;

    if(r->ri_tes.allocated_N == N)
    {    
        r->ri_tes.sim->fPerformSummation(r, r->ri_tes.sim->radau->Qout, r->ri_tes.sim->radau->Pout, 
                                        r->ri_tes.sim->radau->dQ, r->ri_tes.sim->radau->dP, 8);
                    
        double * Q_out = r->ri_tes.sim->radau->Qout;
        double * P_out = r->ri_tes.sim->radau->Pout;
        double * m = r->ri_tes.sim->mass;
        for(uint32_t i=1; i < N; i++) // Change index to zero once all inertial frames are supported.
        {
            r->ri_tes.particles_dh[i].x = Q_out[3*i];
            r->ri_tes.particles_dh[i].y = Q_out[3*i+1];
            r->ri_tes.particles_dh[i].z = Q_out[3*i+2];
            r->ri_tes.particles_dh[i].vx = P_out[3*i]/m[i];    
            r->ri_tes.particles_dh[i].vy = P_out[3*i+1]/m[i];
            r->ri_tes.particles_dh[i].vz = P_out[3*i+2]/m[i];

            r->ri_tes.particles_dh[i].m = m[i];
        }       
        reb_transformations_democraticheliocentric_to_inertial_posvel(r->particles, r->ri_tes.particles_dh, r->N, r->N); 
    }         
}

void reb_integrator_tes_reset(struct reb_simulation* r){
    // Need to think about this properly - should this also go inside the part2 loop? When is reset called?
    // UniversalVars_Free();
    // dhem_Free();
    // Radau_Free();
    // Simulation_Free();    
}
