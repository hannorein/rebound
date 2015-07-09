/**
 * @file 	problem.c
 * @brief 	How to use unique ids to identify particles
 * @author 	Daniel Tamayo <d.tamayo@utoronto.ca>
 * @detail  This example shows how to assign IDs to particles, and demonstrates different 
 * options for removing particles from the simulation.
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Dan Tamayo
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
	// IDs can be set to any integer
	for (int i=0;i<10;i++){
		struct particle p;
		p.x  = i; 		p.y  = i;	 	p.z  = i;
		p.vx = i; 		p.vy = i;	 	p.vz = i;
		p.ax = 0; 			p.ay = 0; 			p.az = 0;
		p.m  = i;
		p.ID = i; // IDs only work if you export PARTICLEIDS=1 in the Makefile
		particles_add(p); 
	}

	printf("Initial IDs:\n");
	print_IDs();

	int success;
	int keepSorted = 0;
	printf("\nTry to remove index 3...\n");
	success = particles_remove(3, keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
	print_IDs();
	printf("Because keepSorted = 0, last particle replaced removed particle and indices got scrambled:\n\n");

	keepSorted = 1;
	printf("Try to remove index 6 while preserving the order with keepSorted=1...\n");
	success = particles_remove(6, keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
	print_IDs();

	printf("\nWe can also remove particles by ID.  Try to remove ID=5...\n");
	success = particles_remove_ID(5, keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
	print_IDs();
	
	printf("\nIf we try to remove an index > N or an ID that doesn't exist, we get a warning and no particle is removed:\n");
	printf("Try to remove index 15...\n");
	success = particles_remove(15, keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
	printf("Try to remove ID=3...\n");
	success = particles_remove_ID(3, keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
	exit(0);
}

void print_IDs(void){
	printf("IDs = ");
	for (int i=0;i<N;i++){
		printf("%d ", particles[i].ID);
	}
	printf("\n");
}

void problem_output(){
}

void problem_finish(){
}
