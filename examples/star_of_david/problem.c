/**
 * @file 	problem.c
 * @brief 	Example problem: Star of David.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example uses the IAS15 integrator
 * to integrate the "Star od David", a four body system consisting of two
 * binaries orbiting each other. Note that the time is running backwards,
 * which illustrates that IAS15 can handle both forward and backward in time
 * integrations. The initial conditions are by Robert Vanderbei. For more 
 * information see http://www.princeton.edu/%7Ervdb/WebGL/New.html
 *
 * @section 	LICENSE
 * Copyright (c) 2014 Hanno Rein, Robert Vanderbei
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
#include "output.h"
#include "tools.h"
#include "particle.h"
#include "boundaries.h"



void problem_init(int argc, char* argv[]){
	init_boxwidth(6.); 			
	dt = -1;

	struct particle p = {.m = 1., .z = 0., .vz = 0.};
	
	p.x =  -1.842389804706855; p.y =  -1.063801316823613; 
	p.vx =  -0.012073765486548; p.vy =   0.021537467220014; 
	particles_add(p);

	p.x =  -0.689515464218133; p.y =  -0.398759403276399; 
	p.vx =   0.637331229856386; p.vy =  -1.103822313621890; 
	particles_add(p);
	
	p.x =   0.689515464218133; p.y =   0.398759403276399; 
	p.vx =  -0.637331229856386; p.vy =   1.103822313621890; 
	particles_add(p);
	
	p.x =   1.842389804706855; p.y =   1.063801316823613; 
	p.vx =   0.012073765486548; p.vy =  -0.021537467220014; 
	particles_add(p);
	
	tools_move_to_center_of_momentum();
}

void problem_output(){
	if (output_check(1.)){
		output_timing();
	}
}

void problem_finish(){
}
