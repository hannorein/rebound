/**
 * @file 	boundaries.c
 * @brief 	Implementation of periodic boundary conditions. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	The code supports different boundary conditions.
 * This file implements simple periodic boundary conditions.
 * 
 * 
 * @section LICENSE
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
#include "particle.h"
#include "integrator.h"
#include "main.h"
#include "boundaries.h"
#include "tree.h"

int nghostx = 1;
int nghosty = 1;
int nghostz = 1;

void boundaries_check(){
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		while(particles[i].x>boxsize_x/2.){
			particles[i].x -= boxsize_x;
		}
		while(particles[i].x<-boxsize_x/2.){
			particles[i].x += boxsize_x;
		}
		while(particles[i].y>boxsize_y/2.){
			particles[i].y -= boxsize_y;
		}
		while(particles[i].y<-boxsize_y/2.){
			particles[i].y += boxsize_y;
		}
		while(particles[i].z>boxsize_z/2.){
			particles[i].z -= boxsize_z;
		}
		while(particles[i].z<-boxsize_z/2.){
			particles[i].z += boxsize_z;
		}
	}
}

struct ghostbox boundaries_get_ghostbox(int i, int j, int k){
	struct ghostbox gb;
	gb.shiftx = boxsize_x*(double)i;
	gb.shifty = boxsize_y*(double)j;
	gb.shiftz = boxsize_z*(double)k;
	gb.shiftvx = 0;
	gb.shiftvy = 0;
	gb.shiftvz = 0;
	return gb;
}


