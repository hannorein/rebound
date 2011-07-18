#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "collisions.h"
#include "collision_resolve.h"
#include "main.h"
#include "tree.h"
#include "boundaries.h"

struct collision* collisions = NULL;
int collisions_NMAX = 0;
int collisions_N = 0;
double nearest_r;
double collisions_max_r;
int nearest_pt;
struct ghostbox nearest_gb;

void collisions_add(int p1, int p2, struct ghostbox gb);
void tree_get_nearest_neighbour_in_cell(int pt, struct ghostbox gb, struct cell* c);

void collisions_search(){
	tree_update();
	int nghostxcol = (nghostx>1?1:nghostx);
	int nghostycol = (nghosty>1?1:nghosty);
	int nghostzcol = (nghostz>1?1:nghostz);
	for (int i=0;i<N;i++){
		nearest_r = boxsize_max;
		nearest_pt = -1;
		for (int gbx=-nghostxcol; gbx<=nghostxcol; gbx++){
		for (int gby=-nghostycol; gby<=nghostycol; gby++){
		for (int gbz=-nghostzcol; gbz<=nghostzcol; gbz++){
			struct ghostbox gb = get_ghostbox(gbx,gby,gbz);
			for (int ri=0;ri<root_nx*root_ny*root_nz;ri++){
				struct cell* rootcell = root[ri];
				if (rootcell!=NULL){
					tree_get_nearest_neighbour_in_cell(i,gb,rootcell);
				}
			}
		}
		}
		}
		if (nearest_pt==-1) continue;
		struct particle p1 = particles[i];
		struct particle p2 = particles[nearest_pt];
		double dx = p1.x - (p2.x+nearest_gb.shiftx); 
		double dy = p1.y - (p2.y+nearest_gb.shifty); 
		double dz = p1.z - (p2.z+nearest_gb.shiftz); 
		double sr = p1.r + p2.r; 
		double r2 = dx*dx+dy*dy+dz*dz;
		if (r2>sr*sr) continue;	// not overlapping
		double dvx = p1.vx - (p2.vx+nearest_gb.shiftvx); 
		double dvy = p1.vy - (p2.vy+nearest_gb.shiftvy); 
		double dvz = p1.vz - (p2.vz+nearest_gb.shiftvz); 
		if (dvx*dx + dvy*dy + dvz*dz >0) continue; // not approaching
		collisions_add(i,nearest_pt, nearest_gb);
	}
}

void collisions_add(int p1, int p2, struct ghostbox gb){
	if (collisions_NMAX<=collisions_N){
		collisions_NMAX += 32;
		collisions = realloc(collisions,sizeof(struct collision)*collisions_NMAX);
	}
	collisions[collisions_N].p1 = p1;
	collisions[collisions_N].p2 = p2;
	collisions[collisions_N].gb = gb;
	collisions_N++;
}

void collisions_resolve(){
	for (int i=0;i<collisions_N;i++){
		collisions_resolve_single(collisions[i]);
	}
	collisions_N=0;
}

void tree_get_nearest_neighbour_in_cell(int pt, struct ghostbox gb, struct cell* c){
	struct particle p1 = particles[pt];
	double dx = p1.x - (c->x+gb.shiftx);
	double dy = p1.y - (c->y+gb.shifty);
	double dz = p1.z - (c->z+gb.shiftz);
	double r2 = dx*dx + dy*dy + dz*dz;
	double crit1 = p1.r + collisions_max_r + c->w*0.86602540;
	if (r2 < crit1*crit1){
		if (c->oct!=NULL){
			for (int i=0;i<8;i++){
				if (c->oct[i]!=NULL){
					tree_get_nearest_neighbour_in_cell(pt,gb,c->oct[i]);
				}
			}
		}
	}
	if (c->pt>=0 && c->pt !=pt){
		struct particle p2 = particles[c->pt];
		double dx = p1.x - (p2.x+gb.shiftx);
		double dy = p1.y - (p2.y+gb.shifty);
		double dz = p1.z - (p2.z+gb.shiftz);
		double r2 = dx*dx+dy*dy+dz*dz;
		if (r2<nearest_r*nearest_r){
			nearest_r = sqrt(r2);
			nearest_pt = c->pt;
			nearest_gb = gb;
		}
	}
}
