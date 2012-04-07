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
 * The sweeping direction in this implementation is phi. This can be used
 * for narrow rings, such as in the example 'spreading_ring'.
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
int	N_collisions		= 0;

static inline double min(double a, double b){ return (a>b)?b:a;}
static inline double max(double a, double b){ return (b>a)?b:a;}
static inline double sgn(const double a){ return (a>=0 ? 1. : -1); }

/** 
 * This function checks if two particles colliding during one drift step.
 * @param pt1 Particle 1. 
 * @param pt2 Particle 2. 
 * @param proci Processor id (OpenMP) for this collision.
 * @param crossing Flag that is one if one of the particles crosses a boundary in this timestep.
 */
void detect_collision_of_pair(int pt1, int pt2, int proci, int crossing);

/**
 * Structure that stores a start or end point of a particle trajectory.
 */
struct phivalue {
	double 	phi;		// position along sweep axis
	int 	inout;		// start or endpoint
	int	nphi;		
	int 	crossing;	// crosses boundary
	int 	pt;		// particle
};

/**
 * Structure that contains a list of xvalues.
 */
struct phivaluelist {
	struct phivalue* phivalues;
	int 	N;		/**< Current array size. */
	int 	Nmax;		/**< Maximum array size before realloc() is needed. */
};
struct  phivaluelist* sweepphi;	/**< Pointers to the SWEEPY list of each processor. */

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
 * Adds a line to the SWEEPY array of processor proci.
 */
void add_line_to_phivsublist(double phi1, double phi2, int pt, int n, int proci, int crossing){
	int N = sweepphi[proci].N;
	
	if (N+2>sweepphi[proci].Nmax){
		sweepphi[proci].Nmax 	+= 1024;
		sweepphi[proci].phivalues	= (struct phivalue*)realloc(sweepphi[proci].phivalues,sweepphi[proci].Nmax*sizeof(struct phivalue));
	}

	sweepphi[proci].phivalues[N].phi 		= phi1;
	sweepphi[proci].phivalues[N].pt 		= pt;
	sweepphi[proci].phivalues[N].nphi 		= n;
	sweepphi[proci].phivalues[N].inout 		= 0;
	sweepphi[proci].phivalues[N].crossing 	= crossing;
	sweepphi[proci].phivalues[N+1].phi 		= phi2;
	sweepphi[proci].phivalues[N+1].pt 		= pt;
	sweepphi[proci].phivalues[N+1].nphi 		= n;
	sweepphi[proci].phivalues[N+1].inout	= 1;
	sweepphi[proci].phivalues[N+1].crossing	= crossing;

	sweepphi[proci].N += 2;
}

/**
 * Adds a line to the SWEEPY array and checks for crossings of processor boundaries.
 */
void add_line_to_phivlist(double phi1, double phi2, int pt, int n, int crossing){
	int prociphi1 = (int)(floor( (phi1/(2.*M_PI)+0.5) *(double)sweeps_proc));// %sweeps.phivlists;
	int prociphi2 = (int)(floor( (phi2/(2.*M_PI)+0.5) *(double)sweeps_proc));// %sweeps.phivlists;
	if (prociphi2>=sweeps_proc){
		prociphi2 = sweeps_proc-1;
	}
	if (prociphi1<0){
		prociphi1 = 0;
	}

	if (prociphi1!=prociphi2){
		double b = -M_PI+2.*M_PI/(double)sweeps_proc*(double)prociphi2; 
		add_line_to_phivsublist(phi1,b,pt,n,prociphi1,1);
		add_line_to_phivsublist(b,phi2,pt,n,prociphi2,1);
	}else{
		add_line_to_phivsublist(phi1,phi2,pt,n,prociphi1,crossing);
	}
}

/**
 * Adds a line to the SWEEPY array and checks for crossings of simulation boundaries.
 */
void add_to_phivlist(double phi1, double phi2, int pt){
	double phimin, phimax;
	if (phi1 < phi2){
		phimin = phi1;
		phimax = phi2;
	}else{
		phimin = phi2;
		phimax = phi1;
	}
	double radius = particles[pt].r*1.0001; //Safety factor to avoid floating point issues.
	phimin -= radius;
	phimax += radius;

	if (phimin<-M_PI){
		add_line_to_phivlist(phimin+2.*M_PI,M_PI,pt,1,1);
		add_line_to_phivlist(-M_PI,phimax,pt,0,1);
		return;
	}
	if (phimax>M_PI){
		add_line_to_phivlist(-M_PI,phimax-2.*M_PI,pt,-1,1);
		add_line_to_phivlist(phimin,M_PI,pt,0,1);
		return;
	}
	add_line_to_phivlist(phimin,phimax,pt,0,0);
}

/**
 * Compares the phi position of two phivalues.
 */
int compare_phivalue (const void * a, const void * b){
	const double diff = ((struct phivalue*)a)->phi - ((struct phivalue*)b)->phi;
	if (diff > 0) return 1;
	if (diff < 0) return -1;
	return 0;
}

/**
 * Compares the phi position of two particles.
 */
int compare_particle (const void * a, const void * b){
	const double diff = atan2(((struct particle*)a)->y,((struct particle*)a)->x) - atan2(((struct particle*)b)->y,((struct particle*)b)->x);
	if (diff > 0) return 1;
	if (diff < 0) return -1;
	return 0;
}

/**
 * Sorts the array phivl with insertion sort.
 */
void collisions_sweep_insertionsort_phivaluelist(struct phivaluelist* phivl){
	struct phivalue* phiv = phivl->phivalues;
	int _N = phivl->N;
	for(int j=1;j<_N;j++){
		struct phivalue key = phiv[j];
		int i = j - 1;
		while(i >= 0 && phiv[i].phi > key.phi){
		    phiv[i+1] = phiv[i];
		    i--;
		}
		phiv[i+1] = key;
	}
}

/**
 * Sorts the particle array with insertion sort.
 */
void collisions_sweep_insertionsort_particles(){
	for(int j=1+N_collisions;j<N;j++){
		struct particle key = particles[j];
		double keyphi = atan2(particles[j].y,particles[j].x);
		int i = j - 1;
		while(i >= N_collisions && atan2(particles[i].y,particles[i].x) > keyphi){
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
		sweeps_proc 	= omp_get_max_threads();
#endif // OPENMP
		sweepphi	= (struct phivaluelist*) calloc(sweeps_proc,sizeof(struct phivaluelist));
		clist		= (struct collisionlist*)calloc(sweeps_proc,sizeof(struct collisionlist));
#ifndef TREE
		// Sort particles according to their phi position to speed up sorting of lines.
		// Initially the particles are not pre-sorted, thus qsort is faster than insertionsort.
		// Note that this rearranges particles and will cause problems if the particle id is used elsewhere.
		qsort (&(particles[N_collisions]), N-N_collisions, sizeof(struct particle), compare_particle);
	}else{
		// Keep particles sorted according to their phi position to speed up sorting of lines.
		collisions_sweep_insertionsort_particles();
#endif //TREE
	}
	for (int i=N_collisions;i<N;i++){
		double phi  = atan2(particles[i].y,particles[i].x);
		if (phi != phi) continue;
		double r = sqrt(particles[i].x*particles[i].x + particles[i].y*particles[i].y);
		double w = (particles[i].x*particles[i].vy - particles[i].y*particles[i].vx) / r;
		if (w != w) continue;
		double oldphi = phi-0.5*dt*w-collisions_max_r/r*2.*M_PI;	
		double newphi = phi+0.5*dt*w+collisions_max_r/r*2.*M_PI;	
		add_to_phivlist(oldphi,newphi,i);
	}
	
#pragma omp parallel for schedule (static,1)
	for (int proci=0;proci<sweeps_proc;proci++){
		struct phivaluelist* sweepphii = &(sweepphi[proci]);
#ifdef TREE
		// Use quicksort when there is a tree. Particles are not pre-sorted.
		qsort (sweepphii->phivalues, sweepphii->N, sizeof(struct phivalue), compare_phivalue);
#else //TREE 
		// Use insertionsort when there is a tree. Particles are pre-sorted.
		collisions_sweep_insertionsort_phivaluelist(sweepphii);	
#endif //TREE
		
		// SWEEPL: List of lines intersecting the plane.
		struct phivaluelist sweepl = {NULL,0,0};

		for (int i=0;i<sweepphii->N;i++){
			struct phivalue phiv = sweepphii->phivalues[i];
			if (phiv.inout == 0){
				// Add event if start of line
				if (sweepl.N>=sweepl.Nmax){
					sweepl.Nmax +=32;
		 			sweepl.phivalues = realloc(sweepl.phivalues,sizeof(struct phivalue)*sweepl.Nmax); 
				}
				sweepl.phivalues[sweepl.N] = phiv;
				// Check for collisions with other particles in SWEEPL
				for (int k=0;k<sweepl.N;k++){
					int p1 = phiv.pt;
					int p2 = sweepl.phivalues[k].pt;
					if (p1==p2) continue;
					int gbnphi = phiv.nphi;
					if (sweepl.phivalues[k].nphi!=0){
						if (sweepl.phivalues[k].nphi==phiv.nphi) continue;
						int tmp = p1;
						p1 = p2;
						p2 = tmp;
						gbnphi = sweepl.phivalues[k].nphi;
					}
					detect_collision_of_pair(p1,p2,proci,sweepl.phivalues[k].crossing||phiv.crossing);
				}
				sweepl.N++;
			}else{
				// Remove event if end of line
				for (int j=0;j<sweepl.N;j++){
					if (sweepl.phivalues[j].pt == phiv.pt){
						sweepl.N--;
						sweepl.phivalues[j] = sweepl.phivalues[sweepl.N];
						j--;
						break;
					}
				}
			}
		}
		free(sweepl.phivalues);
	}

}

void detect_collision_of_pair(int pt1, int pt2, int proci, int crossing){
	struct particle* p1 = &(particles[pt1]);
	struct particle* p2 = &(particles[pt2]);
	double x  = p1->x   - p2->x;
	double y  = p1->y   - p2->y;
	double z  = p1->z   - p2->z;
	double vx = p1->vx  - p2->vx;
	double vy = p1->vy  - p2->vy;
	double vz = p1->vz  - p2->vz;

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
			c1.gb = boundaries_get_ghostbox(0,0,0);
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
		sweepphi[proci].N = 0;
	}
#ifdef OPENMP
	omp_destroy_lock(&boundarylock);
#endif //OPENMP
}

