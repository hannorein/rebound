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


#define n  3
static double Q[n*3];
static double P[n*3];

static SIMULATION * sim;

void reb_integrator_tes_part1(struct reb_simulation* r){} // unused

void reb_integrator_tes_part2(struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    uint32_t N = r->N;

    // for(uint32_t i=0;i<N;i++) 
    // {
    //     Q[3*i] = particles[i].x;
    //     Q[3*i+1] = particles[i].y;
    //     Q[3*i+2] = particles[i].z;

    //     P[3*i] = particles[i].vx;
    //     P[3*i+1] = particles[i].vy;
    //     P[3*i+2] = particles[i].vz;        
    // }

    // remove this when an actual step is performed (enables testing from python for now)
    // for(uint32_t i=0;i<N;i++) 
    // {
    //     particles[i].x = 0;
    //     particles[i].y = 0;
    //     particles[i].z = 0;
    //     particles[i].vx = 0;
    //     particles[i].vy = 0;
    //     particles[i].vz = 0;
    //     particles[i].ax = 0;
    //     particles[i].ay = 0;
    //     particles[i].az = 0;
    // }        

    // update timestep
	r->t+=r->dt;
	r->dt_last_done = r->dt;

    // update the particle arrays - dont need to call this all the time. Does the main rebound routine handle this when required?
    // reb_integrator_tes_synchronize(r);
}

void reb_integrator_tes_synchronize(struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    uint32_t N = r->N;

    // // Do I need to convert away from DH coords here?
    // for(uint32_t i=0; i < N; i++) 
    // {
    //     particles[i].x = Q[3*i];
    //     particles[i].y = Q[3*i+1];
    //     particles[i].z = Q[3*i+2];

    //     // This need to be converted back to velocity
    //     particles[i].vx = P[3*i];
    //     particles[i].vy = P[3*i+1];
    //     particles[i].vz = P[3*i+2];
    // }
}

void reb_integrator_tes_reset(struct reb_simulation* r){
	// Do nothing.
}

// @todo this is temporary code and needs to be removed once unit testing is wrapped around
void reb_integrator_tes_init(struct reb_simulation* r)
{
    struct reb_particle* const particles = r->particles;
    uint32_t N = r->N;

        // Allocate memory for input arrays.
    double * Q = (double *)malloc(sizeof(double)*3*r->N);
    double * V = (double *)malloc(sizeof(double)*3*r->N);
    double * m = (double *)malloc(sizeof(double)*r->N);

    sim = Simulation_Init(N);

    // Default values - these need changing to be proper
    double t0 = 0;
    double period = r->ri_tes.orbital_period;
    double orbits = 100;
    uint32_t output_spacing = 0;
    uint32_t output_samples = r->ri_tes.output_samples;
    double rTol = r->ri_tes.epsilon;
    double aTol = 1;
    double rectisPerOrbit = r->ri_tes.recti_per_orbit;
    double dQcutoff = r->ri_tes.dq_max;
    double dPcutoff = 1;
    double hInitial = r->dt;
    double timeOut = 1e11;
    char * outputFile = "temp_output.txt";  

    for(uint32_t i=0;i<N;i++) 
    {
        Q[3*i] = particles[i].x;
        Q[3*i+1] = particles[i].y;
        Q[3*i+2] = particles[i].z;

        V[3*i] = particles[i].vx;
        V[3*i+1] = particles[i].vy;
        V[3*i+2] = particles[i].vz;        

        m[i] = particles[i].m*r->G; // need to think carefully about what the do with G (probs change to rebound style to make use of helper functions).
    }

    // Perform initialisation of our integrator objects.
    Sim_AddInitialConditions(Q, V, m);

    // Store configuration values in the simulation.
    sim->t0 = t0;
    sim->tEnd = period*orbits;
    sim->hInitial = hInitial;
    sim->orbits = orbits;
    sim->aTol = aTol;
    sim->rTol = rTol;
    sim->rectisPerOrbit = rectisPerOrbit;
    sim->orbits = orbits;
    sim->period = period;
    sim->dQcutoff = dQcutoff;
    sim->dPcutoff = dPcutoff;
    sim->outputFile = outputFile;
    sim->timeOut = timeOut;
    sim->output_spacing = output_spacing;
    sim->output_samples = output_samples;

    UniversalVars_Init(sim);
    dhem_Init(sim, period/rectisPerOrbit, 9);
    Radau_Init(sim);
    // Perform the integration.
    Radau_integrate();

    // Store the output data.
    double * Q_out = sim->radau->Qout;
    double * P_out = sim->radau->Pout;
    for(uint32_t i=0; i < N; i++) // Do I need to convert away from DH coords here?
    {
        particles[i].x = Q_out[3*i];
        particles[i].y = Q_out[3*i+1];
        particles[i].z = Q_out[3*i+2];

        // This need to be converted back to velocity
        particles[i].vx = P_out[3*i];
        particles[i].vy = P_out[3*i+1];
        particles[i].vz = P_out[3*i+2];
    }
    
    // Clean up after onesself.
    // UniversalVars_Free();
    // dhem_Free();
    // Radau_Free();
    // Simulation_Free();
}