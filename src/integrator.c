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
#include "rebound.h"
#include "gravity.h"
#include "output.h"
#include "integrator.h"
#include "integrator_whfast.h"
#include "integrator_ias15.h"
#include "integrator_leapfrog.h"
#include "integrator_sei.h"
#include "integrator_wh.h"
#include "integrator_hybrid.h"

void reb_integrator_part1(struct reb_simulation* r){
	switch(r->integrator){
		case REB_INTEGRATOR_IAS15:
			reb_integrator_ias15_part1(r);
			break;
		case REB_INTEGRATOR_WH:
			reb_integrator_wh_part1(r);
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
		case REB_INTEGRATOR_HYBRID:
			reb_integrator_hybrid_part1(r);
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
		case REB_INTEGRATOR_WH:
			reb_integrator_wh_part2(r);
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
		case REB_INTEGRATOR_HYBRID:
			reb_integrator_hybrid_part2(r);
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
		case REB_INTEGRATOR_WH:
			reb_integrator_wh_synchronize(r);
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
		case REB_INTEGRATOR_HYBRID:
			reb_integrator_hybrid_synchronize(r);
			break;
		default:
			break;
	}
}

void reb_integrator_reset(struct reb_simulation* r){
	r->integrator = REB_INTEGRATOR_IAS15;
	r->gravity_ignore_10 = 0;
	reb_integrator_ias15_reset(r);
	reb_integrator_wh_reset(r);
	reb_integrator_leapfrog_reset(r);
	reb_integrator_sei_reset(r);
	reb_integrator_whfast_reset(r);
	reb_integrator_hybrid_reset(r);
}

void reb_update_acceleration(struct reb_simulation* r){
	// This should probably go elsewhere
	PROFILING_STOP(PROFILING_CAT_INTEGRATOR)
	PROFILING_START()
	reb_calculate_acceleration(r);
	if (r->N_var){
		reb_calculate_acceleration_var(r);
	}
	if (r->additional_forces) r->additional_forces(r);
	PROFILING_STOP(PROFILING_CAT_GRAVITY)
	PROFILING_START()
}

