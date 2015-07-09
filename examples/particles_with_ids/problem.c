/**
 * @file 	problem.c
 * @brief 	Example problem: solar system.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example integrates all planets of the Solar
 * System. The data comes from the NASA HORIZONS system. 
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Hanno Rein, Shangfei Liu, Dave Spiegel
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
#include "integrator.h"
#include "integrator_whfast.h"

void print_IDs(void);

void problem_init(int argc, char* argv[]){
	// Setup constants
	tmax		= 0.;			
	// Just demonstrating IDs, so set initial conditions arbitrarily (and IDs in the order they are added)
	for (int i=0;i<10;i++){
		struct particle p;
		p.x  = i; 		p.y  = i;	 	p.z  = i;
		p.vx = i; 		p.vy = i;	 	p.vz = i;
		p.ax = 0; 			p.ay = 0; 			p.az = 0;
		p.m  = i;
		p.ID = i; // IDs only work if you export PARTICLEIDS = 1 in the Makefile
		particles_add(p); 
	}

	print_IDs();

	int success;
	int index = 3;
	int keepSorted = 0;
	success = particles_remove(index, keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
	print_IDs(); // if keepSorted = 0, last particle replaces removed particle, and indices get scrambled
	
	keepSorted = 1;
	success = particles_remove(7, keepSorted);
	print_IDs();

	success = particles_remove_ID(5, keepSorted);
	print_IDs();

	success = particles_remove(15, keepSorted);
	success = particles_remove_ID(3, keepSorted);
}

void print_IDs(void){
	printf("IDs = ");
	for (int i=0;i<10;i++){
		printf("%d ", particles[i].ID);
	}
	printf("\n");
}

void problem_output(){
}

void problem_finish(){
}
