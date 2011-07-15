#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "collisions.h"
#include "main.h"
#include "tree.h"
#include "boundaries.h"

double collisions_max_r;

struct collision{
	struct particle* p1;
	struct particle* p2;
	struct ghostbox gb;
} collision;

struct collision* collisions = NULL;
int collisions_NMAX = 0;
int collisions_N = 0;
double coefficient_of_restitution = 1;
double minimum_collision_velocity = 0;
void collisions_add(struct particle* p1, struct particle* p2, struct ghostbox gb);
void tree_get_nearest_neighbour_in_cell(struct particle* p, struct ghostbox gb, struct cell* c);

double nearest_r;
struct particle* nearest_p;
struct ghostbox nearest_gb;

void collisions_search(){
	if (root==NULL){
		tree_init();
	} else {
		check_boundaries();
		tree_check_moved_particles();
		tree_update(root);
	}
	int nghostxcol = (nghostx>1?1:nghostx);
	int nghostycol = (nghosty>1?1:nghosty);
	int nghostzcol = (nghostz>1?1:nghostz);
	for (int i=0;i<N;i++){
		struct particle* p1 = &(particles[i]);
		nearest_r = boxsize_max;
		nearest_p = NULL;
		for (int gbx=-nghostxcol; gbx<=nghostxcol; gbx++){
		for (int gby=-nghostycol; gby<=nghostycol; gby++){
		for (int gbz=-nghostzcol; gbz<=nghostzcol; gbz++){
			struct ghostbox gb = get_ghostbox(gbx,gby,gbz);
			for (int ri=0;ri<root_nx;ri++){
			for (int rj=0;rj<root_ny;rj++){
			for (int rk=0;rk<root_nz;rk++){
				if (gbx== 1 && ri!=0) continue;
				if (gbx==-1 && ri!=0) continue;
				if (gby== 1 && rj!=0) continue;
				if (gby==-1 && rj!=0) continue;
				if (gbz== 1 && rk!=0) continue;
				if (gbz==-1 && rk!=0) continue;
				int index = (rk*root_ny+rj)*root_nx+ri;
				tree_get_nearest_neighbour_in_cell(p1,gb,root[index]);
			}
			}
			}
		}
		}
		}
		if (nearest_p==NULL) continue;
		double dx = p1->x - (nearest_p->x+nearest_gb.shiftx); 
		double dy = p1->y - (nearest_p->y+nearest_gb.shifty); 
		double dz = p1->z - (nearest_p->z+nearest_gb.shiftz); 
		double sr = p1->r + nearest_p->r; 
		double r2 = dx*dx+dy*dy+dz*dz;
		if (r2>sr*sr) continue;	// not overlapping
		double dvx = p1->vx - (nearest_p->vx-nearest_gb.shiftvx); 
		double dvy = p1->vy - (nearest_p->vy-nearest_gb.shiftvy); 
		double dvz = p1->vz - (nearest_p->vz-nearest_gb.shiftvz); 
		if (dvx*dx + dvy*dy + dvz*dz >0) continue; // not approaching
		collisions_add(p1,nearest_p, nearest_gb);
	}
}

void collisions_add(struct particle* p1, struct particle* p2, struct ghostbox gb){
	if (collisions_NMAX<=collisions_N){
		collisions_NMAX += 32;
		collisions = realloc(collisions,sizeof(struct collision)*collisions_NMAX);
	}
	collisions[collisions_N].p1 = p1;
	collisions[collisions_N].p2 = p2;
	collisions[collisions_N].gb = gb;
	collisions_N++;
}

void collisions_resolve_single(struct collision c){
	struct particle p1 = *(c.p1);
	struct particle p2 = *(c.p2);
	// The following line would disallow multiple collisons per timestep.
//	if (p1.lastcollision==t || p2.lastcollision==t) return;
	double m21  = p1.m  /  p2.m; 
	double x21  = p1.x  - (p2.x+c.gb.shiftx); 
	double y21  = p1.y  - (p2.y+c.gb.shifty); 
	double z21  = p1.z  - (p2.z+c.gb.shiftz); 
	double vx21 = p1.vx - (p2.vx-c.gb.shiftvx); 
	double vy21 = p1.vy - (p2.vy-c.gb.shiftvy); 
	double vz21 = p1.vz - (p2.vz-c.gb.shiftvz); 
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
	c.p2->vx -=	m21*dvx2n;
	c.p2->vy -=	m21*dvy2nn;
	c.p2->vz -=	m21*dvz2nn;
	c.p2->lastcollision = t;
	c.p1->vx +=	dvx2n; 
	c.p1->vy +=	dvy2nn; 
	c.p1->vz +=	dvz2nn; 
	c.p1->lastcollision = t;

}

void collisions_resolve(){
	for (int i=0;i<collisions_N;i++){
		collisions_resolve_single(collisions[i]);
	}
	collisions_N=0;
}

void tree_get_nearest_neighbour_in_cell(struct particle* p, struct ghostbox gb, struct cell* c){
	double dx = p->x - (c->x+gb.shiftx);
	double dy = p->y - (c->y+gb.shifty);
	double dz = p->z - (c->z+gb.shiftz);
	double r2 = dx*dx+dy*dy+dz*dz;
	if (r2<nearest_r*nearest_r+c->w*c->w*3.+2.*nearest_r+2.*1.732*c->w){
		if (c->oct!=NULL){
			for (int i=0;i<8;i++){
				if (c->oct[i]!=NULL){
					tree_get_nearest_neighbour_in_cell(p,gb,c->oct[i]);
				}
			}
		}
	}
	if (c->pt!=NULL && c->pt!=p){
		double dx = p->x - (c->pt->x+gb.shiftx);
		double dy = p->y - (c->pt->y+gb.shifty);
		double dz = p->z - (c->pt->z+gb.shiftz);
		double r2 = dx*dx+dy*dy+dz*dz;
		if (r2<nearest_r*nearest_r){
			nearest_r = sqrt(r2);
			nearest_p = c->pt;
			nearest_gb = gb;
		}
	}
}
