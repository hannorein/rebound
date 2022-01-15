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

        Simulation_Init(r, N);

        // Convert from inertial to dh coords.
        reb_transformations_inertial_to_democraticheliocentric_posvel(particles, r->ri_tes.particles_dh, r->N, r->N);

        // Store the COM and it's velocity as they are treated differently by TES.
        r->ri_tes.COM[0] = r->ri_tes.particles_dh[0].x;
        r->ri_tes.COM[1] = r->ri_tes.particles_dh[0].y;
        r->ri_tes.COM[2] = r->ri_tes.particles_dh[0].z;

        r->ri_tes.COM_dot[0] = r->ri_tes.particles_dh[0].vx;
        r->ri_tes.COM_dot[1] = r->ri_tes.particles_dh[0].vy;
        r->ri_tes.COM_dot[2] = r->ri_tes.particles_dh[0].vz;        
        
        for(uint32_t i=1;i<N;i++) 
        {
            r->ri_tes.mass[i] =     r->ri_tes.particles_dh[i].m;
            r->ri_tes.Q_dh[3*i] =   r->ri_tes.particles_dh[i].x;
            r->ri_tes.Q_dh[3*i+1] = r->ri_tes.particles_dh[i].y;
            r->ri_tes.Q_dh[3*i+2] = r->ri_tes.particles_dh[i].z;
            r->ri_tes.P_dh[3*i] =   r->ri_tes.particles_dh[i].vx*r->ri_tes.particles_dh[i].m;
            r->ri_tes.P_dh[3*i+1] = r->ri_tes.particles_dh[i].vy*r->ri_tes.particles_dh[i].m;
            r->ri_tes.P_dh[3*i+2] = r->ri_tes.particles_dh[i].vz*r->ri_tes.particles_dh[i].m; 
        }
        r->ri_tes.mass[0] = particles[0].m; // Keep mass[0] as stellar mass.

        UniversalVars_Init(r);
        dhem_Init(r, r->ri_tes.orbital_period/r->ri_tes.recti_per_orbit, 9);
        dhem_InitialiseOsculatingOrbits(r, r->ri_tes.Q_dh, r->ri_tes.P_dh, r->t);
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
        dhem_PerformSummation(r, r->ri_tes.radau->Qout, r->ri_tes.radau->Pout, 
                                        r->ri_tes.radau->dQ, r->ri_tes.radau->dP, 8);
                    
        double * Q_out = r->ri_tes.radau->Qout;
        double * P_out = r->ri_tes.radau->Pout;
        double * m = r->ri_tes.mass;

        for(uint32_t i=1; i < N; i++)
        {
            r->ri_tes.particles_dh[i].x = Q_out[3*i];
            r->ri_tes.particles_dh[i].y = Q_out[3*i+1];
            r->ri_tes.particles_dh[i].z = Q_out[3*i+2];
            r->ri_tes.particles_dh[i].vx = P_out[3*i]/m[i];    
            r->ri_tes.particles_dh[i].vy = P_out[3*i+1]/m[i];
            r->ri_tes.particles_dh[i].vz = P_out[3*i+2]/m[i];

            r->ri_tes.particles_dh[i].m = m[i];
        }       

        // Update the output with the COM including drift.
        r->ri_tes.particles_dh[0].x = r->ri_tes.COM[0];
        r->ri_tes.particles_dh[0].y = r->ri_tes.COM[1];
        r->ri_tes.particles_dh[0].z = r->ri_tes.COM[2];

        r->ri_tes.particles_dh[0].vx = r->ri_tes.COM_dot[0];
        r->ri_tes.particles_dh[0].vy = r->ri_tes.COM_dot[1];
        r->ri_tes.particles_dh[0].vz = r->ri_tes.COM_dot[2];  

        r->ri_tes.particles_dh[0].m = r->ri_tes.rhs->mTotal;

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

void reb_integrator_tes_allocate_memory(struct reb_simulation* r)
{
    Simulation_Init(r, r->N);
    r->ri_tes.particles_dh = (struct particles*)malloc(sizeof(struct reb_particle)*r->N);
    UniversalVars_Init(r);
    dhem_Init(r, r->ri_tes.orbital_period/r->ri_tes.recti_per_orbit, 9);
    Radau_Init(r);          
}
