/**
 * @file    integrator_saba.c
 * @brief   SABA integrator family (Laskar and Robutel 2001).
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details This file implements the family of symplectic integrators
 * of Laskar and Robutel (2001). 
 * 
 * @section LICENSE
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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include "rebound.h"
#include "particle.h"
#include "tools.h"
#include "gravity.h"
#include "integrator.h"
#include "integrator_whfast.h"
#include "integrator_saba.h"

#define MAX(a, b) ((a) < (b) ? (b) : (a))   ///< Returns the maximum of a and b
#define MIN(a, b) ((a) > (b) ? (b) : (a))   ///< Returns the minimum of a and b

// Some coefficients appear multiple times to simplify the loop structures. 
const static double reb_saba_c[4][4] = {
        {0.5, 0., 0., 0.}, // SABA1
        {0.2113248654051871177454256097490212721762, 0.5773502691896257645091487805019574556476, 0., 0.}, // SABA2
        {0.1127016653792583114820734600217600389167, 0.3872983346207416885179265399782399610833,  0.3872983346207416885179265399782399610833, 0.}, // SABA3
        {0.06943184420297371238802675555359524745214, 0.2605776340045981552106403648947824089476, 0.3399810435848562648026657591032446872006, 0.2605776340045981552106403648947824089476}, // SABA4
}; 
const static double reb_saba_d[4][4] = {
        {1., 0., 0., 0.},
        {0.5, 0.5, 0., 0.},
        {0.2777777777777777777777777777777777777778, 0.4444444444444444444444444444444444444444,0.2777777777777777777777777777777777777778, 0.},
        {0.1739274225687269286865319746109997036177, 0.3260725774312730713134680253890002963823, 0.3260725774312730713134680253890002963823, 0.1739274225687269286865319746109997036177},
};
const static double reb_saba_cc[4] = {
        0.08333333333333333333333333333333333333333, // SABA1
        0.01116454968463011276968973577058865137738, // SABA2
        0.005634593363122809402267823769797538671562, // SABA3
        0.003396775048208601331532157783492144, // SABA4
}; 
    

static void reb_saba_corrector_step(struct reb_simulation* r, double cc){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_particle* const p_j = ri_whfast->p_jh;
	struct reb_particle* const particles = r->particles;
    const int N = r->N;
    switch (r->ri_saba.corrector){
        case REB_SABA_CORRECTOR_MODIFIEDKICK: 
            // Calculate normal kick
            reb_transformations_jacobi_to_inertial_pos(particles, p_j, particles, N);
            reb_update_acceleration(r);
            // Calculate jerk
            reb_whfast_calculate_jerk(r);

            for (int i=0; i<N; i++){
                const double prefact = r->dt*r->dt;
                particles[i].ax = prefact*p_j[i].ax; 
                particles[i].ay = prefact*p_j[i].ay; 
                particles[i].az = prefact*p_j[i].az; 
            }
            reb_whfast_interaction_step(r,cc*r->dt);
            break;
        case REB_SABA_CORRECTOR_LAZY: 
            {
            // Need temporary array to store old positions
            if (ri_whfast->allocated_Ntemp != N){
                ri_whfast->allocated_Ntemp = N;
                ri_whfast->p_temp = realloc(ri_whfast->p_temp,sizeof(struct reb_particle)*N);
            }
            struct reb_particle* p_temp = ri_whfast->p_temp;

            // Calculate normal kick
            reb_transformations_jacobi_to_inertial_pos(particles, p_j, particles, N);
            reb_update_acceleration(r);
            reb_transformations_inertial_to_jacobi_acc(particles, p_j, particles, N);

            // make copy of original positions and accelerations
            memcpy(p_temp,p_j,r->N*sizeof(struct reb_particle));

            // WHT96 Eq 10.6
            const double prefac1 = r->dt*r->dt/12.; 
            for (unsigned int i=1;i<N;i++){
                p_j[i].x += prefac1 * p_temp[i].ax;
                p_j[i].y += prefac1 * p_temp[i].ay;
                p_j[i].z += prefac1 * p_temp[i].az;
            }
           
            // recalculate kick 
            reb_transformations_jacobi_to_inertial_pos(particles, p_j, particles, N);
            reb_update_acceleration(r);
            reb_transformations_inertial_to_jacobi_acc(particles, p_j, particles, N);

            const double prefact = cc*r->dt*12.;
            for (unsigned int i=1;i<N;i++){
                // Lazy implementer's commutator
                p_j[i].vx += prefact*(p_j[i].ax - p_temp[i].ax);
                p_j[i].vy += prefact*(p_j[i].ay - p_temp[i].ay);
                p_j[i].vz += prefact*(p_j[i].az - p_temp[i].az);
                // reset positions
                p_j[i].x = p_temp[i].x;
                p_j[i].y = p_temp[i].y;
                p_j[i].z = p_temp[i].z;
            }
            }
            break;
        default:
            return;
    };
}

void reb_integrator_saba_part1(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_simulation_integrator_saba* const ri_saba = &(r->ri_saba);
    const int k = ri_saba->k;
    const int corrector = ri_saba->corrector;
    if (r->var_config_N>0){
        reb_error(r, "Variational particles are not supported in the SABA integrator.");
        return; 
    }
    if (ri_whfast->coordinates!=REB_WHFAST_COORDINATES_JACOBI){
        reb_error(r, "SABA integrator requires ri_whfast.coordinates to be set to Jacobi coordinates.");
        return; 
    }
    if (k>4 || k==0){
        reb_error(r, "SABA is only implemented up to SABA4.");
        return; 
    }
    if (corrector>2){
        reb_error(r, "SABA corrector setting needs to be 0 (off), 1 (modified kick), or 2 (lazy implementer's method).");
        return; 
    }
    if (corrector){
        // Force Jacobi terms to be calculated in reb_update_acceleration if corrector is used
        r->gravity = REB_GRAVITY_JACOBI;
    }else{
        // Otherwise can do either way
        r->gravity_ignore_terms = 1;
    }
    if (reb_integrator_whfast_init(r)){
        // Non recoverable error occured.
        return;
    }
    
    // Only recalculate Jacobi coordinates if needed
    if (ri_saba->safe_mode || ri_whfast->recalculate_coordinates_this_timestep){
        reb_integrator_whfast_from_inertial(r);
        ri_whfast->recalculate_coordinates_this_timestep = 0;
    }
    if (corrector){
        if (ri_saba->is_synchronized){
            reb_saba_corrector_step(r, reb_saba_cc[k-1]);
        }else{
            reb_saba_corrector_step(r, 2.*reb_saba_cc[k-1]);
        }
        // First half DRIFT step
        reb_whfast_kepler_step(r, reb_saba_c[k-1][0]*r->dt);   
        reb_whfast_com_step(r, reb_saba_c[k-1][0]*r->dt);
    }else{
        if (ri_saba->is_synchronized){
            // First half DRIFT step
            reb_whfast_kepler_step(r, reb_saba_c[k-1][0]*r->dt);   
            reb_whfast_com_step(r, reb_saba_c[k-1][0]*r->dt);
        }else{
            // Combined DRIFT step
            reb_whfast_kepler_step(r, 2.*reb_saba_c[k-1][0]*r->dt);   
            reb_whfast_com_step(r, 2.*reb_saba_c[k-1][0]*r->dt);
        }
    }

    reb_integrator_whfast_to_inertial(r);
}

void reb_integrator_saba_synchronize(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_simulation_integrator_saba* const ri_saba = &(r->ri_saba);
    int k = ri_saba->k;
    if (ri_saba->is_synchronized == 0){
        const int N = r->N;
        if (ri_saba->corrector){
            // Drift already done, just need corrector
            reb_saba_corrector_step(r, reb_saba_cc[k-1]);
        }else{
            reb_whfast_kepler_step(r, reb_saba_c[k-1][0]*r->dt);
            reb_whfast_com_step(r, reb_saba_c[k-1][0]*r->dt);
        }
        reb_transformations_jacobi_to_inertial_posvel(r->particles, ri_whfast->p_jh, r->particles, N);
        ri_saba->is_synchronized = 1;
    }
}

void reb_integrator_saba_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_simulation_integrator_saba* const ri_saba = &(r->ri_saba);
    struct reb_particle* restrict const particles = r->particles;
    int k = ri_saba->k;
    const int N = r->N;
    if (ri_whfast->p_jh==NULL){
        // Non recoverable error occured earlier. 
        // Skipping rest of integration to avoid segmentation fault.
        return;
    }
    
    reb_whfast_interaction_step(r, reb_saba_d[k-1][0]*r->dt);
  
    for(int i=1;i<k;i++){
        reb_whfast_kepler_step(r, reb_saba_c[k-1][i]*r->dt);   
        reb_whfast_com_step(r, reb_saba_c[k-1][i]*r->dt);
        reb_transformations_jacobi_to_inertial_pos(particles, ri_whfast->p_jh, particles, N);
        reb_update_acceleration(r);
        reb_whfast_interaction_step(r, reb_saba_d[k-1][i]*r->dt);
    } 

    if (ri_saba->corrector){
        // Always need to do drift step if correctors are turned on
        reb_whfast_kepler_step(r, reb_saba_c[k-1][0]*r->dt);
        reb_whfast_com_step(r, reb_saba_c[k-1][0]*r->dt);
    }

    ri_saba->is_synchronized = 0;
    if (ri_saba->safe_mode){
        reb_integrator_saba_synchronize(r);
    }
    
    r->t+=r->dt;
    r->dt_last_done = r->dt;
}
    
void reb_integrator_saba_reset(struct reb_simulation* const r){
    struct reb_simulation_integrator_saba* const ri_saba = &(r->ri_saba);
    ri_saba->k = 1;
    ri_saba->corrector = REB_SABA_CORRECTOR_NONE;
    ri_saba->safe_mode = 1;
    ri_saba->is_synchronized = 1;
    reb_integrator_whfast_reset(r);
}
