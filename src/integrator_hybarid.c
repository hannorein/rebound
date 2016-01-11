/**
 * @file 	integrator.c
 * @brief 	Leap-frog integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	This file implements the leap-frog integration scheme.  
 * This scheme is second order accurate, symplectic and well suited for 
 * non-rotating coordinate systems. Note that the scheme is formally only
 * first order accurate when velocity dependent forces are present.
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
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
#include "integrator_ias15.h"
#include "integrator_whfast.h"


void reb_integrator_hybarid_additional_forces(struct reb_simulation* mini){
}

void reb_integrator_hybarid_part1(struct reb_simulation* r){
    if (r->ri_hybarid.mini == NULL){
        r->ri_hybarid.mini = reb_create_simulation();
        r->ri_hybarid.mini->integrator = REB_INTEGRATOR_IAS15;
        r->ri_hybarid.mini->additional_forces = reb_integrator_hybarid_additional_forces;
    }
    reb_integrator_whfast_part1(r);
}
void reb_integrator_hybarid_part2(struct reb_simulation* r){
    reb_integrator_whfast_part2(r);

    //
    //
    //
    // Run mini
    //
    // Check for encounters
    //   if new encounters, copy from global to mini
    //
    // Copy positions from mini to global
    //
    // Store the pos to prev_pos
    //
}
	
void reb_integrator_hybarid_synchronize(struct reb_simulation* r){
	// Do nothing.
    reb_integrator_whfast_synchronize(r);
}

void reb_integrator_hybarid_reset(struct reb_simulation* r){
	// Do nothing.
    reb_integrator_whfast_reset(r);
}
