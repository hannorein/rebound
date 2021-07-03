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
#include "simulation.h"
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

        r->ri_tes.sim = Simulation_Init(N);

        // Default values - these need changing to be proper
        double t0 = 0;
        double orbits = 100;

        for(uint32_t i=0;i<N;i++) 
        {
            r->ri_tes.sim->Q0[3*i] = particles[i].x;
            r->ri_tes.sim->Q0[3*i+1] = particles[i].y;
            r->ri_tes.sim->Q0[3*i+2] = particles[i].z;
            r->ri_tes.sim->V0[3*i] = particles[i].vx;
            r->ri_tes.sim->V0[3*i+1] = particles[i].vy;
            r->ri_tes.sim->V0[3*i+2] = particles[i].vz;        
            r->ri_tes.sim->mass[i] = particles[i].m*r->G; // need to think carefully about what the do with G (probs change to rebound style to make use of helper functions).
        }

        // Store configuration values in the simulation. (should be able to remove this section entirely eventually)
        r->ri_tes.sim->t0 = t0;
        r->ri_tes.sim->tEnd = r->ri_tes.orbital_period*r->ri_tes.orbits;
        r->ri_tes.sim->hInitial = r->dt;
        r->ri_tes.sim->aTol = 1;
        r->ri_tes.sim->rTol = r->ri_tes.epsilon;
        r->ri_tes.sim->rectisPerOrbit = r->ri_tes.recti_per_orbit;
        r->ri_tes.sim->period = r->ri_tes.orbital_period;
        r->ri_tes.sim->dQcutoff = r->ri_tes.dq_max;
        r->ri_tes.sim->dPcutoff = 1;
        r->ri_tes.sim->timeOut = 1e11;

        UniversalVars_Init(r->ri_tes.sim);
        dhem_Init(r->ri_tes.sim, r->ri_tes.orbital_period/r->ri_tes.recti_per_orbit, 9);
        Radau_Init(r->ri_tes.sim);  
    }

    double dt_new = Radau_SingleStep(r->t, r->dt, r->dt_last_done);

    // update timestep
	r->t+=r->dt;
	r->dt_last_done = r->dt;
    r->dt = dt_new;
}

void reb_integrator_tes_synchronize(struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    uint32_t N = r->N;

    r->ri_tes.sim->fPerformSummation(r->ri_tes.sim->radau->Qout, r->ri_tes.sim->radau->Pout, 
                                     r->ri_tes.sim->radau->dQ, r->ri_tes.sim->radau->dP, 8);
                 
    double * Q_out = r->ri_tes.sim->radau->Qout;
    double * P_out = r->ri_tes.sim->radau->Pout;
    for(uint32_t i=0; i < N; i++) // Do I need to convert away from DH coords here?
    {
        // Probably OK to just do Qosc+dQ and Posc+dP rather than having a separate function here.
        particles[i].x = Q_out[3*i];
        particles[i].y = Q_out[3*i+1];
        particles[i].z = Q_out[3*i+2];

        // These need to be converted back to velocity
        particles[i].vx = P_out[3*i];
        particles[i].vy = P_out[3*i+1];
        particles[i].vz = P_out[3*i+2];
    }                
}

void reb_integrator_tes_reset(struct reb_simulation* r){
    // Need to think about this properly - should this also go inside the part2 loop?
    // UniversalVars_Free();
    // dhem_Free();
    // Radau_Free();
    // Simulation_Free();    
}
