/**
 * @file 	integrator.c
 * @brief 	Integration schemes.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the leap-frog integration scheme.  
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
#include "main.h"
#include "gravity.h"
#include "problem.h"
#include "output.h"
#include "integrator.h"
#include "integrator_whfast.h"
#include "integrator_ias15.h"
#include "integrator_leapfrog.h"
#include "integrator_sei.h"
#include "integrator_wh.h"
#include "integrator_hybrid.h"

integrator_t integrator = IAS15;
unsigned int integrator_force_is_velocitydependent = 0;

void integrator_part1(void){
	switch(integrator){
		case IAS15:
			integrator_ias15_part1();
			break;
		case WH:
			integrator_wh_part1();
			break;
		case LEAPFROG:
			integrator_leapfrog_part1();
			break;
		case SEI:
			integrator_sei_part1();
			break;
		case WHFAST:
			integrator_whfast_part1();
			break;
		case HYBRID:
			integrator_hybrid_part1();
			break;
		default:
			break;
	}
}

void integrator_part2(void){
	switch(integrator){
		case IAS15:
			integrator_ias15_part2();
			break;
		case WH:
			integrator_wh_part2();
			break;
		case LEAPFROG:
			integrator_leapfrog_part2();
			break;
		case SEI:
			integrator_sei_part2();
			break;
		case WHFAST:
			integrator_whfast_part2();
			break;
		case HYBRID:
			integrator_hybrid_part2();
			break;
		default:
			break;
	}
}
	
void integrator_synchronize(void){
	switch(integrator){
		case IAS15:
			integrator_ias15_synchronize();
			break;
		case WH:
			integrator_wh_synchronize();
			break;
		case LEAPFROG:
			integrator_leapfrog_synchronize();
			break;
		case SEI:
			integrator_sei_synchronize();
			break;
		case WHFAST:
			integrator_whfast_synchronize();
			break;
		case HYBRID:
			integrator_hybrid_synchronize();
			break;
		default:
			break;
	}
}

void integrator_reset(void){
	integrator = IAS15;
	gravity_ignore_10 = 0;
	integrator_force_is_velocitydependent = 0;
	integrator_ias15_reset();
	integrator_wh_reset();
	integrator_leapfrog_reset();
	integrator_sei_reset();
	integrator_whfast_reset();
	integrator_hybrid_reset();
}

void integrator_update_acceleration(void){
	PROFILING_STOP(PROFILING_CAT_INTEGRATOR)
	PROFILING_START()
	gravity_calculate_acceleration();
	if (N_megno){
		gravity_calculate_variational_acceleration();
	}
	if (problem_additional_forces) problem_additional_forces();
	if (problem_additional_forces_with_parameters) problem_additional_forces_with_parameters(particles,t,dt,G,N,N_megno);
	PROFILING_STOP(PROFILING_CAT_GRAVITY)
	PROFILING_START()
}
