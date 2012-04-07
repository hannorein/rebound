/**
 * @file 	collisions.c
 * @brief 	Collision search using a line sweep algorithm, O(N log(N)).
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	The routines in this file implement a collision detection
 * method called line sweep. It is very fast if all dimensions except one 
 * are small. The algorithm is similar to the original algorithm proposed
 * by Bentley & Ottmann (1979) but does not maintain a binary search tree.
 * This is much faster as long as the number of particle trajectories
 * currently intersecting the plane is small.
 *
 * The sweeping direction in this implementation is x.
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
#include "tree.h"
#include "boundaries.h"
#ifdef OPENMP
#include <omp.h>
#endif


double 	collisions_max_r	= 0;
double 	collisions_max2_r	= 0;
int	sweeps_proc		= 1;	/**< Number of processors used for seeping algorithm. */
int 	sweeps_init_done 	= 0;	/**< Used for initialisation of data structures. */

static inline double min(double a, double b){ return (a>b)?b:a;}
static inline double max(double a, double b){ return (b>a)?b:a;}
static inline double sgn(const double a){ return (a>=0 ? 1. : -1); }

/** 
 * This function checks if two particles colliding during one drift step.
 * @param pt1 Particle 1. 
 * @param pt2 Particle 2. 
 * @param proci Processor id (OpenMP) for this collision.
 * @param crossing Flag that is one if one of the particles crosses a boundary in this timestep.
 * @param ghostbox Ghostbox used in this collision.
 */
void detect_collision_of_pair(int pt1, int pt2, int proci, int crossing, struct ghostbox gb);

/**
 * Structure that stores a start or end point of a particle trajectory.
 */
struct xvalue {
	double 	x;		// position along sweep axis
	int 	inout;		// start or endpoint
	int	nx;		
	int 	crossing;	// crosses boundary
	int 	pt;		// particle
};

/**
 * Structure that contains a list of xvalues.
 */
struct xvaluelist {
	struct xvalue* xvalues;
	int 	N;		/**< Current array size. */
	int 	Nmax;		/**< Maximum array size before realloc() is needed. */
};
struct  xvaluelist* sweepx;	/**< Pointers to the SWEEPX list of each processor. */

/**
 * Structure that contains a list of collisions.
 */
struct collisionlist {
	struct collision* collisions;
	int 	N;		/**< Current array size. */
	int 	Nmax;		/**< Maximum array size before realloc() is needed. */
};
struct 	collisionlist* clist;	/**< Pointers to the collisions list of each processor. */

/**
 * Adds a line to the SWEEPX array of processor proci.
 */
void add_line_to_xvsublist(double x1, double x2, int pt, int n, int proci, int crossing){
	int N = sweepx[proci].N;
	
	if (N+2>sweepx[proci].Nmax){
		sweepx[proci].Nmax 	+= 1024;
		sweepx[proci].xvalues	= (struct xvalue*)realloc(sweepx[proci].xvalues,sweepx[proci].Nmax*sizeof(struct xvalue));
	}

	sweepx[proci].xvalues[N].x 		= x1;
	sweepx[proci].xvalues[N].pt 		= pt;
	sweepx[proci].xvalues[N].nx 		= n;
	sweepx[proci].xvalues[N].inout 		= 0;
	sweepx[proci].xvalues[N].crossing 	= crossing;
	sweepx[proci].xvalues[N+1].x 		= x2;
	sweepx[proci].xvalues[N+1].pt 		= pt;
	sweepx[proci].xvalues[N+1].nx 		= n;
	sweepx[proci].xvalues[N+1].inout	= 1;
	sweepx[proci].xvalues[N+1].crossing	= crossing;

	sweepx[proci].N += 2;
}

/**
 * Adds a line to the SWEEPX array and checks for crossings of processor boundaries.
 */
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

/**
 * Adds a line to the SWEEPX array and checks for crossings of simulation boundaries.
 */
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

/**
 * Compares the x position of two xvalues.
 */
int compare_xvalue (const void * a, const void * b){
	const double diff = ((struct xvalue*)a)->x - ((struct xvalue*)b)->x;
	if (diff > 0) return 1;
	if (diff < 0) return -1;
	return 0;
}

/**
 * Compares the x position of two particles.
 */
int compare_particle (const void * a, const void * b){
	const double diff = ((struct particle*)a)->x - ((struct particle*)b)->x;
	if (diff > 0) return 1;
	if (diff < 0) return -1;
	return 0;
}

/**
 * Sorts the array xvl with insertion sort.
 */
void collisions_sweep_insertionsort_xvaluelist(struct xvaluelist* xvl){
	struct xvalue* xv = xvl->xvalues;
	int _N = xvl->N;
	for(int j=1;j<_N;j++){
		struct xvalue key = xv[j];
		int i = j - 1;
		while(i >= 0 && xv[i].x > key.x){
		    xv[i+1] = xv[i];
		    i--;
		}
		xv[i+1] = key;
	}
}

/**
 * Sorts the particle array with insertion sort.
 */
void collisions_sweep_insertionsort_particles(){
	for(int j=1;j<N;j++){
		struct particle key = particles[j];
		int i = j - 1;
		while(i >= 0 && particles[i].x > key.x){
		    particles[i+1] = particles[i];
		    i--;
		}
		particles[i+1] = key;
	}
}



void collisions_search(){
	if (sweeps_init_done!=1){
		sweeps_init_done = 1;
#ifdef OPENMP
		sweeps_proc 		= omp_get_max_threads();
#endif // OPENMP
		sweepx		= (struct xvaluelist*)   calloc(sweeps_proc,sizeof(struct xvaluelist));
		clist		= (struct collisionlist*)calloc(sweeps_proc,sizeof(struct collisionlist));
#ifndef TREE
		// Sort particles according to their x position to speed up sorting of lines.
		// Initially the particles are not pre-sorted, thus qsort is faster than insertionsort.
		// Note that this rearranges particles and will cause problems if the particle id is used elsewhere.
		qsort (particles, N, sizeof(struct particle), compare_particle);
	}else{
		// Keep particles sorted according to their x position to speed up sorting of lines.
		collisions_sweep_insertionsort_particles();
#endif //TREE
	}
	for (int i=0;i<N;i++){
		double oldx = particles[i].x-0.5*dt*particles[i].vx;	
		double newx = particles[i].x+0.5*dt*particles[i].vx;	
		add_to_xvlist(oldx,newx,i);
	}
	
	// Precalculate most comonly used ghostboxes
	const struct ghostbox gb00  = boundaries_get_ghostbox(0,0,0);
	const struct ghostbox gb0p1 = boundaries_get_ghostbox(0,1,0);
	const struct ghostbox gb0m1 = boundaries_get_ghostbox(0,-1,0);

#pragma omp parallel for schedule (static,1)
	for (int proci=0;proci<sweeps_proc;proci++){
		struct xvaluelist* sweepxi = &(sweepx[proci]);
#ifdef TREE
		// Use quicksort when there is a tree. Particles are not pre-sorted.
		qsort (sweepxi->xvalues, sweepxi->N, sizeof(struct xvalue), compare_xvalue);
#else //TREE 
		// Use insertionsort when there is a tree. Particles are pre-sorted.
		collisions_sweep_insertionsort_xvaluelist(sweepxi);	
#endif //TREE
		
		// SWEEPL: List of lines intersecting the plane.
		struct xvaluelist sweepl = {NULL,0,0};

		for (int i=0;i<sweepxi->N;i++){
			struct xvalue xv = sweepxi->xvalues[i];
			if (xv.inout == 0){
				// Add event if start of line
				if (sweepl.N>=sweepl.Nmax){
					sweepl.Nmax +=32;
		 			sweepl.xvalues = realloc(sweepl.xvalues,sizeof(struct xvalue)*sweepl.Nmax); 
				}
				sweepl.xvalues[sweepl.N] = xv;
				// Check for collisions with other particles in SWEEPL
				for (int k=0;k<sweepl.N;k++){
					int p1 = xv.pt;
					int p2 = sweepl.xvalues[k].pt;
					if (p1==p2) continue;
					int gbnx = xv.nx;
					if (sweepl.xvalues[k].nx!=0){
						if (sweepl.xvalues[k].nx==xv.nx) continue;
						int tmp = p1;
						p1 = p2;
						p2 = tmp;
						gbnx = sweepl.xvalues[k].nx;
					}
					if (gbnx==0){
						// Use cached ghostboxes if possible
						detect_collision_of_pair(p1,p2,proci,sweepl.xvalues[k].crossing||xv.crossing,gb00);
						detect_collision_of_pair(p1,p2,proci,sweepl.xvalues[k].crossing||xv.crossing,gb0p1);
						detect_collision_of_pair(p1,p2,proci,sweepl.xvalues[k].crossing||xv.crossing,gb0m1);
					}else{
						for (int gbny = -1; gbny<=1; gbny++){
							struct ghostbox gb = boundaries_get_ghostbox(gbnx,gbny,0);
							detect_collision_of_pair(p1,p2,proci,sweepl.xvalues[k].crossing||xv.crossing,gb);
						}
					}
				}
				sweepl.N++;
			}else{
				// Remove event if end of line
				for (int j=0;j<sweepl.N;j++){
					if (sweepl.xvalues[j].pt == xv.pt){
						sweepl.N--;
						sweepl.xvalues[j] = sweepl.xvalues[sweepl.N];
						j--;
						break;
					}
				}
			}
		}
		free(sweepl.xvalues);
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
			time2=time1;
			time1=tmp;
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
			if ( (time1>-dt/2. && time1<dt/2.)) { 
				c->time 	= time1;
			}else{
				c->time 	= 0;
			}
				
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
			particles[c1.p1].x -= c1.time*particles[c1.p1].vx; 
			particles[c1.p1].y -= c1.time*particles[c1.p1].vy; 
			particles[c1.p1].z -= c1.time*particles[c1.p1].vz; 
			particles[c1.p2].x -= c1.time*particles[c1.p2].vx; 
			particles[c1.p2].y -= c1.time*particles[c1.p2].vy; 
			particles[c1.p2].z -= c1.time*particles[c1.p2].vz; 
#ifdef OPENMP
			if (c1.crossing){
				omp_set_lock(&boundarylock);
			}
#endif //OPENMP
			collision_resolve(c1);
#ifdef OPENMP
			if (c1.crossing){
				omp_unset_lock(&boundarylock);
			}
#endif //OPENMP
			particles[c1.p1].x += c1.time*particles[c1.p1].vx; 
			particles[c1.p1].y += c1.time*particles[c1.p1].vy; 
			particles[c1.p1].z += c1.time*particles[c1.p1].vz; 
			particles[c1.p2].x += c1.time*particles[c1.p2].vx; 
			particles[c1.p2].y += c1.time*particles[c1.p2].vy; 
			particles[c1.p2].z += c1.time*particles[c1.p2].vz; 
		}
		clist[proci].N = 0;
		sweepx[proci].N = 0;
	}
#ifdef OPENMP
	omp_destroy_lock(&boundarylock);
#endif //OPENMP
}

