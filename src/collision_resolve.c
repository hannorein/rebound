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

double coefficient_of_restitution = 1;
double minimum_collision_velocity = 0;
double collisions_constant_coefficient_of_restitution_for_velocity(double v);
double (*coefficient_of_restitution_for_velocity) (double) = collisions_constant_coefficient_of_restitution_for_velocity;
double collisions_plog =0;

void collision_resolve_single(struct collision c){
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
#warning TODO: Make sure this ratio is the right way round.
	double m21  = p1.m  /  p2.m; 
	double x21  = p1.x + gb.shiftx  - p2.x; 
	double y21  = p1.y + gb.shifty  - p2.y; 
	double z21  = p1.z + gb.shiftz  - p2.z; 
	double rp   = p1.r+p2.r;
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
	double dvx2 = -(1.0+eps)*vx21nn/(1.0+m21) ;
	if (dvx2<minimum_collision_velocity){
		dvx2 = minimum_collision_velocity;
	}

	// Now we are rotating backwards
	double dvx2n = cphi * dvx2;		
	double dvy2n = sphi * dvx2;		
	double dvy2nn = ctheta * dvy2n;	
	double dvz2nn = stheta * dvy2n;	

	// Log y-momentum change
	collisions_plog += fabs(dvy2nn*p1.m*x21);

	// Applying the changes to the particles.
#ifdef MPI
	if (isloc==1){
#endif // MPI
		particles[c.p2].vx -=	m21*dvx2n;
		particles[c.p2].vy -=	m21*dvy2nn;
		particles[c.p2].vz -=	m21*dvz2nn;
		particles[c.p2].lastcollision = t;
#ifdef MPI
	}
#endif // MPI
	particles[c.p1].vx +=	dvx2n; 
	particles[c.p1].vy +=	dvy2nn; 
	particles[c.p1].vz +=	dvz2nn; 
	particles[c.p1].lastcollision = t;
#endif // COLLISIONS_NONE
}

double collisions_constant_coefficient_of_restitution_for_velocity(double v){
	return coefficient_of_restitution;
}
