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

#ifdef MPI
#error COLLISIONS_DIRECT not compatible with MPI
#endif

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
		for (int i=0;i<N;i++){
			struct particle p1 = particles[i];
			struct ghostbox gb = get_ghostbox(gbx,gby,gbz);
			// Precalculate shift.
			gb.shiftx += p1.x;
			gb.shifty += p1.y;
			gb.shiftz += p1.z;
			gb.shiftvx += p1.vx;
			gb.shiftvy += p1.vy;
			gb.shiftvz += p1.vz;
			for (int j=0;j<N;j++){
				if (i==j) continue;
				struct particle p2 = particles[j];
				double dx = gb.shiftx - p2.x; 
				double dy = gb.shifty - p2.y; 
				double dz = gb.shiftz - p2.z; 
				double sr = p1.r + p2.r; 
				double r2 = dx*dx+dy*dy+dz*dz;
				if (r2>sr*sr) continue;	// not overlapping
				double dvx = gb.shiftvx - p2.vx; 
				double dvy = gb.shiftvy - p2.vy; 
				double dvz = gb.shiftvz - p2.vz; 
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
