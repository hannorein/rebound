/**
 * @file 	integrator.c
 * @brief 	Integration schemes.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	This file implements the leap-frog integration scheme.  
 * This scheme is second order accurate, symplectic and well suited for 
 * non-rotating coordinate systems. Note that the scheme is formally only
 * first order accurate when velocity dependent forces are present.
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Hanno Rein
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
#include "rebound.h"
#include "gravity.h"
#include "output.h"
#include "integrator.h"
#include "integrator_whfast.h"
#include "integrator_saba.h"
#include "integrator_ias15.h"
#include "integrator_mercurius.h"
#include "integrator_leapfrog.h"
#include "integrator_sei.h"
#include "integrator_janus.h"
#include "integrator_eos.h"

void reb_integrator_part1(struct reb_simulation* r){
	switch(r->integrator){
		case REB_INTEGRATOR_IAS15:
			reb_integrator_ias15_part1(r);
			break;
		case REB_INTEGRATOR_LEAPFROG:
			reb_integrator_leapfrog_part1(r);
			break;
		case REB_INTEGRATOR_SEI:
			reb_integrator_sei_part1(r);
			break;
		case REB_INTEGRATOR_WHFAST:
			reb_integrator_whfast_part1(r);
			break;
		case REB_INTEGRATOR_SABA:
			reb_integrator_saba_part1(r);
			break;
		case REB_INTEGRATOR_MERCURIUS:
			reb_integrator_mercurius_part1(r);
			break;
		case REB_INTEGRATOR_JANUS:
			reb_integrator_janus_part1(r);
			break;
		case REB_INTEGRATOR_EOS:
			reb_integrator_eos_part1(r);
			break;
		default:
			break;
	}
}

void reb_integrator_part2(struct reb_simulation* r){
	switch(r->integrator){
		case REB_INTEGRATOR_IAS15:
			reb_integrator_ias15_part2(r);
			break;
		case REB_INTEGRATOR_LEAPFROG:
			reb_integrator_leapfrog_part2(r);
			break;
		case REB_INTEGRATOR_SEI:
			reb_integrator_sei_part2(r);
			break;
		case REB_INTEGRATOR_WHFAST:
			reb_integrator_whfast_part2(r);
			break;
		case REB_INTEGRATOR_SABA:
			reb_integrator_saba_part2(r);
			break;
		case REB_INTEGRATOR_MERCURIUS:
			reb_integrator_mercurius_part2(r);
			break;
		case REB_INTEGRATOR_JANUS:
			reb_integrator_janus_part2(r);
			break;
		case REB_INTEGRATOR_EOS:
			reb_integrator_eos_part2(r);
			break;
        case REB_INTEGRATOR_NONE:
            r->t += r->dt;
            r->dt_last_done = r->dt;
            break;
		default:
			break;
	}
}
	
void reb_integrator_synchronize(struct reb_simulation* r){
	switch(r->integrator){
		case REB_INTEGRATOR_IAS15:
			reb_integrator_ias15_synchronize(r);
			break;
		case REB_INTEGRATOR_LEAPFROG:
			reb_integrator_leapfrog_synchronize(r);
			break;
		case REB_INTEGRATOR_SEI:
			reb_integrator_sei_synchronize(r);
			break;
		case REB_INTEGRATOR_WHFAST:
			reb_integrator_whfast_synchronize(r);
			break;
		case REB_INTEGRATOR_SABA:
			reb_integrator_saba_synchronize(r);
			break;
		case REB_INTEGRATOR_MERCURIUS:
			reb_integrator_mercurius_synchronize(r);
			break;
		case REB_INTEGRATOR_JANUS:
			reb_integrator_janus_synchronize(r);
			break;
		case REB_INTEGRATOR_EOS:
			reb_integrator_eos_synchronize(r);
			break;
		default:
			break;
	}
}

void reb_integrator_init(struct reb_simulation* r){
	switch(r->integrator){
		case REB_INTEGRATOR_SEI:
			reb_integrator_sei_init(r);
			break;
		default:
			break;
	}
}

void reb_integrator_reset(struct reb_simulation* r){
	r->integrator = REB_INTEGRATOR_IAS15;
	r->gravity_ignore_terms = 0;
	reb_integrator_ias15_reset(r);
	reb_integrator_mercurius_reset(r);
	reb_integrator_leapfrog_reset(r);
	reb_integrator_sei_reset(r);
	reb_integrator_whfast_reset(r);
	reb_integrator_saba_reset(r);
	reb_integrator_janus_reset(r);
	reb_integrator_eos_reset(r);
}

void reb_update_acceleration(struct reb_simulation* r){
	// This should probably go elsewhere
	PROFILING_STOP(PROFILING_CAT_INTEGRATOR)
	PROFILING_START()
	reb_calculate_acceleration(r);
	if (r->N_var){
		reb_calculate_acceleration_var(r);
	}
	if (r->additional_forces  && (r->integrator != REB_INTEGRATOR_MERCURIUS || r->ri_mercurius.mode==0)){
        // For Mercurius:
        // Additional forces are only calculated in the kick step, not during close encounter
        if (r->integrator==REB_INTEGRATOR_MERCURIUS){
            // shift pos and velocity so that external forces are calculated in inertial frame
            // Note: Copying avoids degrading floating point performance
            if(r->N>r->ri_mercurius.allocatedN_additionalforces){
                r->ri_mercurius.particles_backup_additionalforces = realloc(r->ri_mercurius.particles_backup_additionalforces, r->N*sizeof(struct reb_particle));
                r->ri_mercurius.allocatedN_additionalforces = r->N;
            }
            memcpy(r->ri_mercurius.particles_backup_additionalforces,r->particles,r->N*sizeof(struct reb_particle)); 
            reb_integrator_mercurius_dh_to_inertial(r);
        }
        r->additional_forces(r);
        if (r->integrator==REB_INTEGRATOR_MERCURIUS){
            struct reb_particle* restrict const particles = r->particles;
            struct reb_particle* restrict const backup = r->ri_mercurius.particles_backup_additionalforces;
            for (int i=0;i<r->N;i++){
                particles[i].x = backup[i].x;
                particles[i].y = backup[i].y;
                particles[i].z = backup[i].z;
                particles[i].vx = backup[i].vx;
                particles[i].vy = backup[i].vy;
                particles[i].vz = backup[i].vz;
            }
        }
    }
	PROFILING_STOP(PROFILING_CAT_GRAVITY)
	PROFILING_START()
}

