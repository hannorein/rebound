/**
 * @file 	problem.c
 * @brief 	Template file for new problems.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
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
#include <sys/time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "communication_mpi.h"
#include "tree.h"
#include "tools.h"

void problem_init(int argc, char* argv[]){
	if (argc>1){ 						// Try to read boxsize from command line
		boxsize = atof(argv[1]);
	}else{
		boxsize = 100;
	}
	init_box();
	
	struct particle pt;
	pt.x 	= 0; 	pt.y 	= 0; 	pt.z 	= 0;
	pt.vx 	= 0; 	pt.vy 	= 0; 	pt.vz 	= 0;
	pt.ax	= 0; 	pt.ay 	= 0; 	pt.az 	= 0;
	pt.m 	= 0;
	particles_add(pt);
}

void problem_inloop(){
}

void problem_output(){
	output_timing();
}

void problem_finish(){
}
