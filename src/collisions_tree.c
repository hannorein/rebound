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
#include "communication_mpi.h"

struct 	collision* collisions 	= NULL;
int 	collisions_NMAX 	= 0;
int 	collisions_N		= 0;
double 	nearest_r2;
double 	collisions_max_r;
int 	nearest_ri;
double 	p1_r;
struct  collision collision_nearest;

void tree_get_nearest_neighbour_in_cell(struct ghostbox gb, int ri, struct cell* c);

void collisions_search(){
	int nghostxcol = (nghostx>1?1:nghostx);
	int nghostycol = (nghosty>1?1:nghosty);
	int nghostzcol = (nghostz>1?1:nghostz);
	for (int i=0;i<N;i++){
		struct particle p1 = particles[i];
		collision_nearest.p1 = i;
		collision_nearest.p2 = -1;
		p1_r = p1.r;
		nearest_r2 = boxsize_max*boxsize_max/4.;
		for (int gbx=-nghostxcol; gbx<=nghostxcol; gbx++){
		for (int gby=-nghostycol; gby<=nghostycol; gby++){
		for (int gbz=-nghostzcol; gbz<=nghostzcol; gbz++){
			struct ghostbox gb = get_ghostbox(gbx,gby,gbz);
			gb.shiftx += p1.x; 
			gb.shifty += p1.y; 
			gb.shiftz += p1.z; 
			gb.shiftvx += p1.vx; 
			gb.shiftvy += p1.vy; 
			gb.shiftvz += p1.vz; 
			for (int ri=0;ri<root_n;ri++){
				struct cell* rootcell = tree_root[ri];
				if (rootcell!=NULL){
					tree_get_nearest_neighbour_in_cell(gb,ri,rootcell);
				}
			}
		}
		}
		}
		if (collision_nearest.p2==-1) continue;
		if (collisions_NMAX<=collisions_N){
			collisions_NMAX += 32;
			collisions = realloc(collisions,sizeof(struct collision)*collisions_NMAX);
		}
		collisions[collisions_N] = collision_nearest;
		collisions_N++;
	}
}

void collisions_resolve(){
	for (int i=0;i<collisions_N;i++){
		collisions_resolve_single(collisions[i]);
	}
	collisions_N=0;
}

void tree_get_nearest_neighbour_in_cell(struct ghostbox gb, int ri, struct cell* c){
	if (c->pt>=0){ 	
		// Leaf node
		int condition 	= 1;
#ifdef MPI
		int isloc	= 1 ;
		isloc = communication_mpi_rootbox_is_local(ri);
		if (isloc==1){
#endif // MPI
			condition = (c->pt != collision_nearest.p1);
#ifdef MPI
		}
#endif // MPI
		if (condition){
			struct particle p2;
#ifdef MPI
			if (isloc==1){
#endif // MPI
				p2 = particles[c->pt];
#ifdef MPI
			}else{
				int root_n_per_node = root_n/mpi_num;
				int proc_id = ri/root_n_per_node;
				p2 = particles_recv[proc_id][c->pt];
			}
#endif // MPI

			double dx = gb.shiftx - p2.x;
			double dy = gb.shifty - p2.y;
			double dz = gb.shiftz - p2.z;
			double r2 = dx*dx+dy*dy+dz*dz;
			if (r2 > nearest_r2) return;
			double dvx = gb.shiftvx - p2.vx;
			double dvy = gb.shiftvy - p2.vy;
			double dvz = gb.shiftvz - p2.vz;
			if (dvx*dx + dvy*dy + dvz*dz >0) return; // not approaching
			nearest_r2 = r2;
			collision_nearest.ri = ri;
			collision_nearest.p2 = c->pt;
			collision_nearest.gb = gb;
		}
	}else{		
		// Not a leaf node
		double dx = gb.shiftx - c->x;
		double dy = gb.shifty - c->y;
		double dz = gb.shiftz - c->z;
		double r2 = dx*dx + dy*dy + dz*dz;
		double rp  = p1_r + collisions_max_r + 0.86602540378443*c->w;
		if (r2 < rp*rp ){
			for (int o=0;o<8;o++){
				struct cell* d = c->oct[o];
				if (d!=NULL){
					tree_get_nearest_neighbour_in_cell(gb,ri,d);
				}
			}
		}
	}
}
