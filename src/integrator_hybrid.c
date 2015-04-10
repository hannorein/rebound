/**
 * @file 	integrator.c
 * @brief 	Hybrid symplectic/IAS15 integration scheme.
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
#include "particle.h"
#include "main.h"
#include "gravity.h"
#include "boundaries.h"
#include "integrator.h"
#include "integrator_mikkola.h"
#include "integrator_ias15.h"

double integrator_hybrid_switch_radius = 10.; // In units of Hill Radii
static double initial_dt = 0;


static double get_min_distance(){
	const double M0 = particles[0].m;
	const int _N_active = ((N_active==-1)?N:N_active)- N_megno;
	const int _N_real   = N - N_megno;
	double min_distance = 1e308;
	for (int i=1; i<_N_active; i++){
	for (int j=i+1; j<_N_real; j++){
		double Msum = particles[i].m+particles[j].m;
		if (Msum == 0.) continue; // Two testparticles
		const double dxi = particles[0].x - particles[i].x;
		const double dyi = particles[0].y - particles[i].y;
		const double dzi = particles[0].z - particles[i].z;
		const double dvxi = particles[0].vx - particles[i].vx;
		const double dvyi = particles[0].vy - particles[i].vy;
		const double dvzi = particles[0].vz - particles[i].vz;
		const double ri = sqrt(dxi*dxi + dyi*dyi + dzi*dzi);
		const double vi2 = (dvxi*dvxi + dvyi*dvyi + dvzi*dvzi);
		const double dxj = particles[0].x - particles[j].x;
		const double dyj = particles[0].y - particles[j].y;
		const double dzj = particles[0].z - particles[j].z;
		const double dvxj = particles[0].vx - particles[j].vx;
		const double dvyj = particles[0].vy - particles[j].vy;
		const double dvzj = particles[0].vz - particles[j].vz;
		const double rj = sqrt(dxj*dxj + dyj*dyj + dzj*dzj);
		const double vj2 = (dvxj*dvxj + dvyj*dvyj + dvzj*dvzj);

		const double dx = particles[i].x - particles[j].x;
		const double dy = particles[i].y - particles[j].y;
		const double dz = particles[i].z - particles[j].z;
		const double r2 = dx*dx + dy*dy + dz*dz;
		const double r = sqrt(r2);
			
		const double mui = G*(M0+particles[i].m);
		const double ai = -mui/( vi2 - 2.*mui/ri );
		const double muj = G*(M0+particles[j].m);
		const double aj = -muj/( vj2 - 2.*muj/rj );
		// Mutual Hill Radius
		const double rHill = pow((particles[i].m+particles[j].m)/(3.*M0),1./3.)*(ai+aj)/2.;
		const double distance = r/rHill;
		if (distance<min_distance){
			min_distance = distance;
		}
	}
	}
	return min_distance;
}


static double distance;
unsigned int integrator_hybrid_mode = 0; // 0 = symplectic; 1 = IAS15
void integrator_hybrid_part1(){
	distance = get_min_distance();
	if (initial_dt==0.){
		initial_dt = dt;
	}
	if (distance<integrator_hybrid_switch_radius){
		if (integrator_hybrid_mode==0){
			integrator_ias15_reset();
		}
		integrator_hybrid_mode = 1;
	}else{
		if (integrator_hybrid_mode==1){
			integrator_mikkola_reset();
			dt = initial_dt;
		}
		integrator_hybrid_mode = 0;
	}
	switch(integrator_hybrid_mode){
		case 0:
			integrator_mikkola_part1();
			break;
		case 1:
			integrator_ias15_part1();
			break;
	}
}
void integrator_hybrid_part2(){
	switch(integrator_hybrid_mode){
		case 0:
			integrator_mikkola_part2();
			break;
		case 1:
			integrator_ias15_part2();
			break;
	}
}
	
void integrator_hybrid_synchronize(){
	//integrator_mikkola_synchronize();
}

void integrator_hybrid_reset(){
	integrator_hybrid_mode = 0;
	integrator_mikkola_reset();
	integrator_ias15_reset();
	integrator_hybrid_switch_radius = 10.;
	initial_dt = 0.;
}
