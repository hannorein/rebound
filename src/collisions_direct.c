#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "collisions.h"
#include "main.h"

// Dummy. No collision search

struct collision{
	int p1;
	int p2;
} collision;

struct collision* collisions = NULL;
int collisions_NMAX = 0;
int collisions_N = 0;
void collisions_add(int i, int j);

void collisions_search(){
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			if (i==j) continue;
			struct particle p1 = particles[i];
			struct particle p2 = particles[j];
			double dx = p1.x - p2.x; 
			double dy = p1.y - p2.y; 
			double dz = p1.z - p2.z; 
			double sr = p1.r + p2.r; 
			double r2 = dx*dx+dy*dy+dz*dz;
			if (r2>sr*sr) continue;	// not overlapping
			double dvx = p1.vx - p2.vx; 
			double dvy = p1.vy - p2.vy; 
			double dvz = p1.vz - p2.vz; 
			if (dvx*dx + dvy*dy + dvz*dz >0) continue; // not approaching
			collisions_add(i,j);
		}
	}
}

void collisions_add(int i, int j){
	if (collisions_NMAX<=collisions_N){
		collisions_NMAX += 32;
		collisions = realloc(collisions,sizeof(struct collision)*collisions_NMAX);
	}
	collisions[collisions_N].p1 = i;
	collisions[collisions_N].p2 = j;
	collisions_N++;
}

void collisions_resolve_single(struct collision c){
	struct particle p1 = particles[c.p1];
	struct particle p2 = particles[c.p2];
	if (p1.lastcollision==t || p2.lastcollision==t) return;
	double m21  = p2.m / p1.m; 
	double x21  = p2.x - p1.x; 
	double y21  = p2.y - p1.y; 
	double z21  = p2.z - p1.z; 
	double vx21 = p2.vx - p1.vx; 
	double vy21 = p2.vy - p1.vy; 
	double vz21 = p2.vz - p1.vz; 
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
	double eps=1.;
	double dvx2 = -(0.5+0.5*eps)*2.0*vx21nn/(1.0+m21) ;

	// Now we are rotating backwards
	double dvx2n = cphi * dvx2;		
	double dvy2n = sphi * dvx2;		
	double dvy2nn = ctheta * dvy2n;	
	double dvz2nn = stheta * dvy2n;	

	// Applying the changes to the particles.
	particles[c.p1].vx -=	m21*dvx2n;
	particles[c.p1].vy -=	m21*dvy2nn;
	particles[c.p1].vz -=	m21*dvz2nn;
	particles[c.p1].lastcollision = t;
	particles[c.p2].vx +=	dvx2n; 
	particles[c.p2].vy +=	dvy2nn; 
	particles[c.p2].vz +=	dvz2nn; 
	particles[c.p2].lastcollision = t;

}

void collisions_resolve(){
	for (int i=0;i<collisions_N;i++){
		collisions_resolve_single(collisions[i]);
	}
	collisions_N=0;
}
