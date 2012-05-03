/**
 * @file 	problem.c
 * @brief 	Example problem: viewer.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This viewer can display data in the form 
 * x, y, z, r. 
 * 
 * @section 	LICENSE
 * Copyright (c) 2012 Hanno Rein, Shangfei Liu
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
#include <sys/time.h>
#include "main.h"
#include "particle.h"
#include "tools.h"

extern int display_pause_sim;
void particles_add_local(struct particle pt);

void problem_init(int argc, char* argv[]){
	dt 	= 0;
	display_pause_sim = 1;
	boxsize	= 100;
	FILE* f;
	if (argc<=1){					
		f = stdin;
		printf("\n\nReading from stdin.\n");
	}else{
		f = fopen(argv[1],"r");
		printf("\n\nReading from file.\n");
	}
	if (f==NULL){
		printf("\n\nERROR. Cannot open file.\n");
		exit(-1);
	}
	// Read in data.
	double x, y, z, r;
	if (fscanf(f,"%lf %lf %lf %lf", &x, &y, &z, &r)==EOF){
		printf("\n\nERROR. Empty file.\n");
		exit(-1);
	}
	double max_x = x, min_x = x;
	double max_y = y, min_y = y;
	double max_z = z, min_z = z;
	init_boxwidth(1);
	while(fscanf(f, "%lf %lf %lf %lf", &x, &y, &z, &r)!=EOF){
		struct particle p = {.x=x, .y=y, .z=z, .r=r};
		particles_add_local(p);
		if (x>max_x) max_x = x;	
		if (y>max_y) max_y = y;	
		if (z>max_z) max_z = z;	
		if (x<min_x) min_x = x;	
		if (y<min_y) min_y = y;	
		if (z<min_z) min_z = z;	
	}
	printf("\nFound %d particles in file.\n",N);
	fclose(f);

	// Resize box.
	boxsize = max_x - min_x;
	if (max_y - min_y > boxsize) boxsize = max_y - min_y;
	if (max_z - min_z > boxsize) boxsize = max_z - min_z;
	boxsize*=1.01;
	boxsize_max = boxsize_x = boxsize_y = boxsize_z = boxsize;
	// Move particles to center.
	for (int i=0;i<N;i++){
		particles[i].x -= (max_x+min_x)/2.;
		particles[i].y -= (max_y+min_y)/2.;
		particles[i].z -= (max_z+min_z)/2.;
	}
}

void problem_inloop(){
}

void problem_output(){
}

void problem_finish(){
}
