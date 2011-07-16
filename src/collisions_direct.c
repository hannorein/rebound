#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "collisions.h"
#include "main.h"
#include "boundaries.h"


// Dummy. No collision search
double collisions_max_r;

struct collision{
	int p1;
	int p2;
	struct ghostbox gb;
} collision;

struct collision* collisions = NULL;
int collisions_NMAX = 0;
int collisions_N = 0;
double coefficient_of_restitution = 1;
double minimum_collision_velocity = 0;
void collisions_add(int i, int j, struct ghostbox gb);

void collisions_search(){
	int nghostxcol = (nghostx>1?1:nghostx);
	int nghostycol = (nghosty>1?1:nghosty);
	int nghostzcol = (nghostz>1?1:nghostz);
	for (int gbx=-nghostxcol; gbx<=nghostxcol; gbx++){
	for (int gby=-nghostycol; gby<=nghostycol; gby++){
	for (int gbz=-nghostzcol; gbz<=nghostzcol; gbz++){
		struct ghostbox gb = get_ghostbox(gbx,gby,gbz);
		for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			if (i==j) continue;
			struct particle p1 = particles[i];
			struct particle p2 = particles[j];
			double dx = p1.x - (p2.x+gb.shiftx); 
			double dy = p1.y - (p2.y+gb.shifty); 
			double dz = p1.z - (p2.z+gb.shiftz); 
			double sr = p1.r + p2.r; 
			double r2 = dx*dx+dy*dy+dz*dz;
			if (r2>sr*sr) continue;	// not overlapping
			double dvx = p1.vx - (p2.vx+gb.shiftvx); 
			double dvy = p1.vy - (p2.vy+gb.shiftvy); 
			double dvz = p1.vz - (p2.vz+gb.shiftvz); 
			if (dvx*dx + dvy*dy + dvz*dz >0) continue; // not approaching
			collisions_add(i,j, gb);
		}
		}
	}
	}
	}
}

void collisions_add(int i, int j, struct ghostbox gb){
	if (collisions_NMAX<=collisions_N){
		collisions_NMAX += 32;
		collisions = realloc(collisions,sizeof(struct collision)*collisions_NMAX);
	}
	collisions[collisions_N].p1 = i;
	collisions[collisions_N].p2 = j;
	collisions[collisions_N].gb = gb;
	collisions_N++;
}

void collisions_resolve_single(struct collision c){
	struct particle p1 = particles[c.p1];
	struct particle p2 = particles[c.p2];
	if (p1.lastcollision==t || p2.lastcollision==t) return;
	double m21  = p1.m  /  p2.m; 
	double x21  = p1.x  - (p2.x+c.gb.shiftx); 
	double y21  = p1.y  - (p2.y+c.gb.shifty); 
	double z21  = p1.z  - (p2.z+c.gb.shiftz); 
	double vx21 = p1.vx - (p2.vx+c.gb.shiftvx); 
	double vy21 = p1.vy - (p2.vy+c.gb.shiftvy); 
	double vz21 = p1.vz - (p2.vz+c.gb.shiftvz); 
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
	double eps= coefficient_of_restitution;
	double dvx2 = -(0.5+0.5*eps)*2.0*vx21nn/(1.0+m21) ;

	if (dvx2<minimum_collision_velocity){
		dvx2 = minimum_collision_velocity;
	}

	// Now we are rotating backwards
	double dvx2n = cphi * dvx2;		
	double dvy2n = sphi * dvx2;		
	double dvy2nn = ctheta * dvy2n;	
	double dvz2nn = stheta * dvy2n;	

	// Applying the changes to the particles.
	particles[c.p2].vx -=	m21*dvx2n;
	particles[c.p2].vy -=	m21*dvy2nn;
	particles[c.p2].vz -=	m21*dvz2nn;
	particles[c.p2].lastcollision = t;
	particles[c.p1].vx +=	dvx2n; 
	particles[c.p1].vy +=	dvy2nn; 
	particles[c.p1].vz +=	dvz2nn; 
	particles[c.p1].lastcollision = t;

}

void collisions_resolve(){
	for (int i=0;i<collisions_N;i++){
		collisions_resolve_single(collisions[i]);
	}
	collisions_N=0;
}
