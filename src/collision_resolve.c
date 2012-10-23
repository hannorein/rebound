/**
 * @file 	collision_resolve.c
 * @brief 	Resolve a single collision. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
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
#include "collision_resolve.h"
#include "main.h"
#include "boundaries.h"
#include "communication_mpi.h"

#define MAX(a, b) ((a) < (b) ? (b) : (a))


double coefficient_of_restitution = 1;
double minimum_collision_velocity = 0;
double collisions_constant_coefficient_of_restitution_for_velocity(double v);
double (*coefficient_of_restitution_for_velocity) (double) = collisions_constant_coefficient_of_restitution_for_velocity;
double 	collisions_plog =0;	/**< Keep track of momentum exchange (used to calculate collisional viscosity in ring systems. */
long	collisions_Nlog =0;	/**< Keep track of Number of collisions. */
void (*collision_resolve) (struct collision) = collision_resolve_hardsphere;

void collision_resolve_hardsphere(struct collision c){
#ifndef COLLISIONS_NONE
	struct particle p1 = particles[c.p1];
	struct particle p2;
#ifdef MPI
	int isloc = communication_mpi_rootbox_is_local(c.ri);
	if (isloc==1){
#endif // MPI
		p2 = particles[c.p2];
#ifdef MPI
	}else{
		int root_n_per_node = root_n/mpi_num;
		int proc_id = c.ri/root_n_per_node;
		p2 = particles_recv[proc_id][c.p2];
	}
#endif // MPI
//	if (p1.lastcollision==t || p2.lastcollision==t) return;
	struct ghostbox gb = c.gb;
	double x21  = p1.x + gb.shiftx  - p2.x; 
	double y21  = p1.y + gb.shifty  - p2.y; 
	double z21  = p1.z + gb.shiftz  - p2.z; 
	double rp   = p1.r+p2.r;
	double oldvyouter;
	if (x21>0){
	 	oldvyouter = p1.vy;
	}else{
		oldvyouter = p2.vy;
	}
	if (rp*rp < x21*x21 + y21*y21 + z21*z21) return;
	double vx21 = p1.vx + gb.shiftvx - p2.vx; 
	double vy21 = p1.vy + gb.shiftvy - p2.vy; 
	double vz21 = p1.vz + gb.shiftvz - p2.vz; 
	if (vx21*x21 + vy21*y21 + vz21*z21 >0) return; // not approaching
	// Bring the to balls in the xy plane.
	// NOTE: this could probabely be an atan (which is faster than atan2)
	double theta = atan2(z21,y21);
	double stheta = sin(theta);
	double ctheta = cos(theta);
	double vy21n = ctheta * vy21 + stheta * vz21;	
	double y21n = ctheta * y21 + stheta * z21;	
	
	// Bring the two balls onto the positive x axis.
	double phi = atan2(y21n,x21);
	double cphi = cos(phi);
	double sphi = sin(phi);
	double vx21nn = cphi * vx21  + sphi * vy21n;		

	// Coefficient of restitution
	double eps= coefficient_of_restitution_for_velocity(vx21nn);
	double dvx2 = -(1.0+eps)*vx21nn;
	double minr = (p1.r>p2.r)?p2.r:p1.r;
	double maxr = (p1.r<p2.r)?p2.r:p1.r;
	double mindv= minr*minimum_collision_velocity;
	double r = sqrt(x21*x21 + y21*y21 + z21*z21);
	mindv *= 1.-(r - maxr)/minr;
	if (mindv>maxr*minimum_collision_velocity)mindv = maxr*minimum_collision_velocity;
	if (dvx2<mindv) dvx2 = mindv;
	// Now we are rotating backwards
	double dvx2n = cphi * dvx2;		
	double dvy2n = sphi * dvx2;		
	double dvy2nn = ctheta * dvy2n;	
	double dvz2nn = stheta * dvy2n;	


	// Applying the changes to the particles.
#ifdef MPI
	if (isloc==1){
#endif // MPI
	const double p2pf = p1.m/(p1.m+p2.m);
	particles[c.p2].vx -=	p2pf*dvx2n;
	particles[c.p2].vy -=	p2pf*dvy2nn;
	particles[c.p2].vz -=	p2pf*dvz2nn;
	particles[c.p2].lastcollision = t;
#ifdef MPI
	}
#endif // MPI
	const double p1pf = p2.m/(p1.m+p2.m);
	particles[c.p1].vx +=	p1pf*dvx2n; 
	particles[c.p1].vy +=	p1pf*dvy2nn; 
	particles[c.p1].vz +=	p1pf*dvz2nn; 
	particles[c.p1].lastcollision = t;
		
	// Return y-momentum change
	if (x21>0){
		collisions_plog += -fabs(x21)*(oldvyouter-particles[c.p1].vy) * p1.m;
		collisions_Nlog ++;
	}else{
		collisions_plog += -fabs(x21)*(oldvyouter-particles[c.p2].vy) * p2.m;
		collisions_Nlog ++;
	}

#endif // COLLISIONS_NONE
}

double collisions_constant_coefficient_of_restitution_for_velocity(double v){
	return coefficient_of_restitution;
}
