#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "collisions.h"
#include "collision_resolve.h"
#include "main.h"
#include "boundaries.h"


double collisions_max_r;

struct collision* collisions = NULL;
int collisions_NMAX = 0;
int collisions_N = 0;
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

void collisions_resolve(){
	for (int i=0;i<collisions_N;i++){
		collisions_resolve_single(collisions[i]);
	}
	collisions_N=0;
}
