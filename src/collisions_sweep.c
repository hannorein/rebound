/**
 * @file 	collisions.c
 * @brief 	Collision search using a line sweep algorithm, O(N log(N)).
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	The routines in this file implement a collision detection
 * method called line sweep. It is very fast for a low number of effective 
 * dimensions (less then three).
 * 
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
#include "collisions.h"
#include "collision_resolve.h"
#include "main.h"
#include "boundaries.h"


// Dummy. No collision search
double 	collisions_max_r	= 0;
int	sweeps_proc;
int 	sweeps_init_done 	= 0;

static inline double min(double a, double b){ return (a>b)?b:a;}
static inline double max(double a, double b){ return (b>a)?b:a;}
static inline double sgn(const double a){ return (a>=0 ? 1. : -1); }

void collisions_resolve_single(struct collision c);
void detect_collision_of_pair(const int pt1, const int pt2, struct ghostbox const gb, int proc, int crossing);


struct xvalue {
	double 	x;		// position along sweep axis
	int 	inout;		// start or endpoint
	int	nx;		
	int 	crossing;	// crosses boundary
	int 	pt;		// particle
};

struct xvaluelist {
	struct xvalue* xvalues;
	int 	N;
	int	_N; 		// array size
};
struct  xvaluelist* restrict xvlist;

struct collisionlist {
	struct collision* collisions;
	int N;
	int _N;
};
struct 	collisionlist* restrict clist;

void init_sweep(){
#ifdef _OPENMP
	sweeps_proc 		= omp_get_max_threads();
#else
	sweeps_proc 		= 1;
#endif
	printf("Optimizing sweep lists for %d processors.\n",sweeps_proc);
	xvlist		= (struct xvaluelist*)calloc(sweeps_proc,sizeof(struct xvaluelist));
	clist		= (struct collisionlist*)calloc(sweeps_proc,sizeof(struct collisionlist));
	for (int i=0;i<sweeps_proc;i++){
		xvlist[i].N		= 0;
		xvlist[i]._N 		= 512;
		xvlist[i].xvalues 	= (struct xvalue*)malloc(xvlist[i]._N*sizeof(struct xvalue));
	}
}

void add_line_to_xvsublist(double x1, double x2, int pt, int n, int i, int crossing){
	int N = xvlist[i].N;
	
	if (N+2>xvlist[i]._N){
		xvlist[i]._N 		+= 1024;
		xvlist[i].xvalues	= (struct xvalue*)realloc(xvlist[i].xvalues,xvlist[i]._N*sizeof(struct xvalue));
		printf("xvlist size now %d\n",xvlist[i]._N);
	}

	xvlist[i].xvalues[N].x 	= x1;
	xvlist[i].xvalues[N].pt 	= pt;
	xvlist[i].xvalues[N].nx 	= n;
	xvlist[i].xvalues[N].inout 	= 0;
	xvlist[i].xvalues[N].crossing 	= crossing;
	xvlist[i].xvalues[N+1].x 	= x2;
	xvlist[i].xvalues[N+1].pt 	= pt;
	xvlist[i].xvalues[N+1].nx 	= n;
	xvlist[i].xvalues[N+1].inout	= 1;
	xvlist[i].xvalues[N+1].crossing= crossing;

	xvlist[i].N += 2;
}

void add_line_to_xvlist(double x1, double x2, int pt, int n, int crossing){
	int ix1 = (int)(floor( (x1/boxsize_x+0.5) *(double)sweeps_proc));// %sweeps.xvlists;
	int ix2 = (int)(floor( (x2/boxsize_x+0.5) *(double)sweeps_proc));// %sweeps.xvlists;
	if (ix2>=sweeps_proc){
		ix2 = sweeps_proc-1;
	}
	if (ix1<0){
		ix1 = 0;
	}

	if (ix1!=ix2){
		double b = -boxsize_x/2.+boxsize_x/(double)sweeps_proc*(double)ix2; 
		add_line_to_xvsublist(x1,b,pt,n,ix1,1);
		add_line_to_xvsublist(b,x2,pt,n,ix2,1);
	}else{
		add_line_to_xvsublist(x1,x2,pt,n,ix1,crossing);
	}
}

void add_to_xvlist(double x1, double x2, int pt){
	double xmin, xmax;
	if (x1 < x2){
		xmin = x1;
		xmax = x2;
	}else{
		xmin = x2;
		xmax = x1;
	}
	const double radius = particles[pt].r;
	xmin -= radius;
	xmax += radius;

	if (xmin<-boxsize_x/2.){
		add_line_to_xvlist(xmin+boxsize_x,boxsize_x/2.,pt,1,1);
		add_line_to_xvlist(-boxsize_x/2.,xmax,pt,0,1);
		return;
	}
	if (xmax>boxsize_x/2.){
		add_line_to_xvlist(-boxsize_x/2.,xmax-boxsize_x,pt,-1,1);
		add_line_to_xvlist(xmin,boxsize_x/2.,pt,0,1);
		return;
	}
	add_line_to_xvlist(xmin,xmax,pt,0,0);
}

int compare_xparticle (const void * a, const void * b){
	const struct particle* x1 = *(struct particle**)a;
	const struct particle* x2 = *(struct particle**)b;
	const double diff = x1->x - x2->x;
	if (diff > 0) return 1;
	if (diff < 0) return -1;
	return 0;
}
int compare_xvalue (const void * a, const void * b){
	const double diff = ((struct xvalue*)a)->x - ((struct xvalue*)b)->x;
	if (diff > 0) return 1;
	if (diff < 0) return -1;
	return 0;
}

void collisions_search(){
	if (sweeps_init_done!=1){
		sweeps_init_done = 1;
		init_sweep();	
	}
	for (int i=0;i<N;i++){
		particles[i].ovx = particles[i].vx;
		particles[i].ovy = particles[i].vy;
		particles[i].ovz = particles[i].vz;
	
		double oldx = particles[i].x-0.5*dt*particles[i].vx;	
		double newx = particles[i].x+0.5*dt*particles[i].vx;	
		add_to_xvlist(oldx,newx,i);
	}
//#pragma omp parallel for
	for (int proci=0;proci<sweeps_proc;proci++){
		struct xvaluelist xvlisti = xvlist[proci];
		qsort (xvlisti.xvalues, xvlisti.N, sizeof(struct xvalue), compare_xvalue);

		struct xvalue** sweeps 	= malloc(sizeof(struct xvalue*)*1024); // Active list. Make the magic number go away.
		int sweeps_N			= 0;

		for (int i=0;i<xvlisti.N;i++){
			struct xvalue* const xv = &(xvlisti.xvalues[i]);
			if (xv->inout == 0){
				// Add event if start of line
				sweeps[sweeps_N] = xv;
				for (int k=0;k<sweeps_N;k++){
					int p1 = xv->pt;
					int p2 = sweeps[k]->pt;
					int gbnx = xv->nx;
					if (sweeps[k]->nx!=0){
						if (sweeps[k]->nx==xv->nx) continue;
						int tmp = p1;
						p1 = p2;
						p2 = tmp;
						gbnx = sweeps[k]->nx;
					}
					for (int gbny = -1; gbny<=1; gbny++){
						struct ghostbox gb = boundaries_get_ghostbox(gbnx,gbny,0);
						detect_collision_of_pair(p1,p2,gb,proci,sweeps[k]->crossing||xv->crossing);
					}
				}
				sweeps_N++;
			}else{
				// Remove event if end of line
				for (int j=0;j<sweeps_N;j++){
					if (sweeps[j]->pt == xv->pt){
						sweeps_N--;
						sweeps[j] = sweeps[sweeps_N];
						j--;
						break;
					}
				}
			}
		}
		free(sweeps);
	}

}

void detect_collision_of_pair(const int pt1, const int pt2, struct ghostbox const gb, int proci, int crossing){
	struct particle* p1 = &(particles[pt1]);
	struct particle* p2 = &(particles[pt2]);
	const double y 	= p1->y	+ gb.shifty	- p2->y;
	const double vy = p1->vy+ gb.shiftvy 	- p2->ovy;

	const double x = p1->x + gb.shiftx	- p2->x;
	const double z = p1->z + gb.shiftz	- p2->z;
	const double vx = p1->ovx + gb.shiftvx	- p2->ovx;
	const double vz = p1->ovz + gb.shiftvz	- p2->ovz;

	const double a = vx*vx + vy*vy + vz*vz;
	const double b = 2.*(vx*x + vy*y + vz*z);
	const double rr = p1->r + p2->r;
	const double c = -rr*rr + x*x + y*y + z*z;

	const double root = b*b-4.*a*c;
	if (root>=0.){
		// Floating point optimized solution of a quadratic equation. Avoids cancelations.
		const double q = -0.5*(b+sgn(b)*sqrt(root));
		const double time2 = q/a;
		const double time1 = c/q;
		//if ( (time1>-dt && time1<0.) || (time2>-dt && time2<0.) || (time1<-dt && time2>0) ){
			// THIS NEEDS IMPROVEMENT #TODO
		if ( (time1>-dt && time1<0.) || (time1<-dt && time2>0) ){
		//printf("Collisiontime: \n\t%f\n\t%f\n",time1,time2);
			double timesave = -dt;
			if (time1>-dt || time2<0){
				if (time1>-dt){
					timesave = time1;
				}else{
					timesave = time2;
				}
			}
			struct collisionlist* const clisti = &(clist[proci]);
			if (clisti->N>=clisti->_N){
				if (clisti->_N==0){
					clisti->_N 		= 1024;
					clisti->collisions	= (struct collision*)malloc(clisti->_N*sizeof(struct collision));
				}else{
					clisti->_N	 	+= 1024;
					clisti->collisions	= (struct collision*)realloc(clisti->collisions,clisti->_N*sizeof(struct collision));
				}
				printf("Collisions array size now: %d\n", clisti->_N);
			}
			struct collision* const c = &(clist->collisions[clisti->N]);
			c->p1		= pt1;
			c->p2		= pt2;
			c->gb	 	= gb;
			c->time 	= timesave;
			c->crossing 	= crossing;
			clisti->N++;
		}
	}
}

void collisions_resolve(){
#ifdef _OPENMP
	omp_lock_t boundarylock;
	omp_init_lock(&boundarylock);
#endif //_OPENMP

//#pragma omp parallel for
	for (int proci=0;proci<sweeps_proc;proci++){
		struct collision* c = clist[proci].collisions;
		int N = clist[proci].N;
	
		// Randomize array.	
		for(int i=0; i<N; i++){
			int j = rand()%N;
			struct collision ctemp = c[i];
			c[i]=c[j];
			c[j]=ctemp;
		}


		for(int i=0; i<N; i++){
			struct collision c1= c[i];
#ifdef _OPENMP
			if (c1.crossing){
				omp_set_lock(&boundarylock);
			}
#endif //_OPENMP
			collision_resolve_single(c1);
#ifdef _OPENMP
			if (c1.crossing){
				omp_unset_lock(&boundarylock);
			}
#endif //_OPENMP
		}
		clist[proci].N = 0;
		xvlist[proci].N = 0;
	}
#ifdef _OPENMP
	omp_destroy_lock(&boundarylock);
#endif //_OPENMP
}

