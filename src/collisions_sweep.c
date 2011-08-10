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
#ifdef OPENMP
#include <omp.h>
#endif


double 	collisions_max_r	= 0;
int	sweeps_proc;
int 	sweeps_init_done 	= 0;

static inline double min(double a, double b){ return (a>b)?b:a;}
static inline double max(double a, double b){ return (b>a)?b:a;}
static inline double sgn(const double a){ return (a>=0 ? 1. : -1); }

void collisions_resolve_single(struct collision c);
void detect_collision_of_pair(int pt1, int pt2, int proci, int crossing, struct ghostbox gb);


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
	int	Nmax; 		// array size
};
struct  xvaluelist* sweepx;

struct collisionlist {
	struct collision* collisions;
	int 	N;
	int 	Nmax;
};
struct 	collisionlist* clist;

void add_line_to_xvsublist(double x1, double x2, int pt, int n, int proci, int crossing){
	int N = sweepx[proci].N;
	
	if (N+2>sweepx[proci].Nmax){
		sweepx[proci].Nmax 		+= 1024;
		sweepx[proci].xvalues	= (struct xvalue*)realloc(sweepx[proci].xvalues,sweepx[proci].Nmax*sizeof(struct xvalue));
	}

	sweepx[proci].xvalues[N].x 		= x1;
	sweepx[proci].xvalues[N].pt 	= pt;
	sweepx[proci].xvalues[N].nx 	= n;
	sweepx[proci].xvalues[N].inout 	= 0;
	sweepx[proci].xvalues[N].crossing 	= crossing;
	sweepx[proci].xvalues[N+1].x 	= x2;
	sweepx[proci].xvalues[N+1].pt 	= pt;
	sweepx[proci].xvalues[N+1].nx 	= n;
	sweepx[proci].xvalues[N+1].inout	= 1;
	sweepx[proci].xvalues[N+1].crossing	= crossing;

	sweepx[proci].N += 2;
}

void add_line_to_xvlist(double x1, double x2, int pt, int n, int crossing){
	int procix1 = (int)(floor( (x1/boxsize_x+0.5) *(double)sweeps_proc));// %sweeps.xvlists;
	int procix2 = (int)(floor( (x2/boxsize_x+0.5) *(double)sweeps_proc));// %sweeps.xvlists;
	if (procix2>=sweeps_proc){
		procix2 = sweeps_proc-1;
	}
	if (procix1<0){
		procix1 = 0;
	}

	if (procix1!=procix2){
		double b = -boxsize_x/2.+boxsize_x/(double)sweeps_proc*(double)procix2; 
		add_line_to_xvsublist(x1,b,pt,n,procix1,1);
		add_line_to_xvsublist(b,x2,pt,n,procix2,1);
	}else{
		add_line_to_xvsublist(x1,x2,pt,n,procix1,crossing);
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
	double radius = particles[pt].r*1.0001; //Safety factor to avoid floating point issues.
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

int compare_xvalue (const void * a, const void * b){
	const double diff = ((struct xvalue*)a)->x - ((struct xvalue*)b)->x;
	if (diff > 0) return 1;
	if (diff < 0) return -1;
	return 0;
}

void collisions_search(){
	if (sweeps_init_done!=1){
		sweeps_init_done = 1;
#ifdef OPENMP
		sweeps_proc 		= omp_get_max_threads();
#else // OPENMP
		sweeps_proc 		= 1;
#endif // OPENMP
		sweepx		= (struct xvaluelist*)   calloc(sweeps_proc,sizeof(struct xvaluelist));
		clist		= (struct collisionlist*)calloc(sweeps_proc,sizeof(struct collisionlist));
	}
	for (int i=0;i<N;i++){
		double oldx = particles[i].x-0.5*dt*particles[i].vx;	
		double newx = particles[i].x+0.5*dt*particles[i].vx;	
		add_to_xvlist(oldx,newx,i);
	}
#pragma omp parallel for schedule (static,1)
	for (int proci=0;proci<sweeps_proc;proci++){
		struct xvaluelist* sweepxi = &(sweepx[proci]);
		qsort (sweepxi->xvalues, sweepxi->N, sizeof(struct xvalue), compare_xvalue);
		
		// SWEEPL list.
		int 		sweepl_N 	= 0;
		int 		sweepl_Nmax 	= 0;
		struct xvalue** sweepl 		= NULL;

		for (int i=0;i<sweepxi->N;i++){
			struct xvalue* xv = &(sweepxi->xvalues[i]);
			if (xv->inout == 0){
				// Add event if start of line
				if (sweepl_N>=sweepl_Nmax){
					sweepl_Nmax +=32;
		 			sweepl= realloc(sweepl,sizeof(struct xvalue*)*sweepl_Nmax); 
				}
				sweepl[sweepl_N] = xv;
				for (int k=0;k<sweepl_N;k++){
					int p1 = xv->pt;
					int p2 = sweepl[k]->pt;
					if (p1==p2) continue;
					int gbnx = xv->nx;
					if (sweepl[k]->nx!=0){
						if (sweepl[k]->nx==xv->nx) continue;
						int tmp = p1;
						p1 = p2;
						p2 = tmp;
						gbnx = sweepl[k]->nx;
					}
					for (int gbny = -1; gbny<=1; gbny++){
						struct ghostbox gb = boundaries_get_ghostbox(gbnx,gbny,0);
						detect_collision_of_pair(p1,p2,proci,sweepl[k]->crossing||xv->crossing,gb);
					}
				}
				sweepl_N++;
			}else{
				// Remove event if end of line
				for (int j=0;j<sweepl_N;j++){
					if (sweepl[j]->pt == xv->pt){
						sweepl_N--;
						sweepl[j] = sweepl[sweepl_N];
						j--;
						break;
					}
				}
			}
		}
		free(sweepl);
	}

}

void detect_collision_of_pair(int pt1, int pt2, int proci, int crossing, struct ghostbox gb){
	struct particle* p1 = &(particles[pt1]);
	struct particle* p2 = &(particles[pt2]);
	double x  = p1->x  + gb.shiftx	- p2->x;
	double y  = p1->y  + gb.shifty	- p2->y;
	double z  = p1->z  + gb.shiftz	- p2->z;
	double vx = p1->vx + gb.shiftvx	- p2->vx;
	double vy = p1->vy + gb.shiftvy - p2->vy;
	double vz = p1->vz + gb.shiftvz	- p2->vz;

	double a = vx*vx + vy*vy + vz*vz;
	double b = 2.*(vx*x + vy*y + vz*z);
	double rr = p1->r + p2->r;
	double c = -rr*rr + x*x + y*y + z*z;

	double root = b*b-4.*a*c;
	if (root>=0.){
		// Floating point optimized solution of a quadratic equation. Avoids cancelations.
		double q = -0.5*(b+sgn(b)*sqrt(root));
		double time1 = c/q;
		double time2 = q/a;
		if (time1>time2){
			double tmp = time2;
			time1=tmp;
			time2=time1;
		}
		if ( (time1>-dt/2. && time1<dt/2.) || (time1<-dt/2. && time2>dt/2.) ){
			struct collisionlist* clisti = &(clist[proci]);
			if (clisti->N>=clisti->Nmax){
				clisti->Nmax	 	+= 1024;
				clisti->collisions	= (struct collision*)realloc(clisti->collisions,clisti->Nmax*sizeof(struct collision));
			}
			struct collision* c = &(clisti->collisions[clisti->N]);
			c->p1		= pt1;
			c->p2		= pt2;
			c->gb	 	= gb;
			c->time 	= time1;
			c->crossing 	= crossing;
			clisti->N++;
		}
	}
}

void collisions_resolve(){
#ifdef OPENMP
	omp_lock_t boundarylock;
	omp_init_lock(&boundarylock);
#endif //OPENMP

#pragma omp parallel for schedule (static,1)
	for (int proci=0;proci<sweeps_proc;proci++){
		struct collision* c = clist[proci].collisions;
		int colN = clist[proci].N;
	
		// Randomize array.	
		for(int i=0; i<colN; i++){
			int j = rand()%colN;
			struct collision ctemp = c[i];
			c[i]=c[j];
			c[j]=ctemp;
		}


		for(int i=0; i<colN; i++){
			struct collision c1= c[i];
			double time = c1.time;
			particles[c1.p1].x += time*particles[c1.p1].vx; 
			particles[c1.p1].y += time*particles[c1.p1].vy; 
			particles[c1.p1].z += time*particles[c1.p1].vz; 
			particles[c1.p2].x += time*particles[c1.p2].vx; 
			particles[c1.p2].y += time*particles[c1.p2].vy; 
			particles[c1.p2].z += time*particles[c1.p2].vz; 
#ifdef OPENMP
			if (c1.crossing){
				omp_set_lock(&boundarylock);
			}
#endif //OPENMP
			collision_resolve_single(c1);
#ifdef OPENMP
			if (c1.crossing){
				omp_unset_lock(&boundarylock);
			}
#endif //OPENMP
			particles[c1.p1].x -= time*particles[c1.p1].vx; 
			particles[c1.p1].y -= time*particles[c1.p1].vy; 
			particles[c1.p1].z -= time*particles[c1.p1].vz; 
			particles[c1.p2].x -= time*particles[c1.p2].vx; 
			particles[c1.p2].y -= time*particles[c1.p2].vy; 
			particles[c1.p2].z -= time*particles[c1.p2].vz; 
		}
		clist[proci].N = 0;
		sweepx[proci].N = 0;
	}
#ifdef OPENMP
	omp_destroy_lock(&boundarylock);
#endif //OPENMP
}

