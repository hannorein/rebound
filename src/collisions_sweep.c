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
	double time;
	int crossing;
} collision;

double coefficient_of_restitution = 1;
double minimum_collision_velocity = 0;

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

struct collisionlist {
	struct collision* collisions;
	int N;
	int _N;
};

struct sweep {
	double	widthx;
	double	widthy;
	struct  xvaluelist* restrict xvlist;
	struct 	collisionlist* restrict clist;
	int	proc;
};

struct sweep sweeps;
int sweeps_init_done = 0;

void init_sweep(){
	sweeps.widthx = boxsize_x;
	sweeps.widthy = boxsize_y;
#ifdef _OPENMP
	sweeps.proc 		= omp_get_max_threads();
#else
	sweeps.proc 		= 1;
#endif
	printf("Optimizing sweep lists for %d processors.\n",sweeps.proc);
	sweeps.xvlist		= (struct xvaluelist*)calloc(sweeps.proc,sizeof(struct xvaluelist));
	sweeps.clist		= (struct collisionlist*)calloc(sweeps.proc,sizeof(struct collisionlist));
	for (int i=0;i<sweeps.proc;i++){
		sweeps.xvlist[i].N		= 0;
		sweeps.xvlist[i]._N 		= 512;
		sweeps.xvlist[i].xvalues 	= (struct xvalue*)malloc(sweeps.xvlist[i]._N*sizeof(struct xvalue));
	}
}

void add_line_to_xvsublist(double x1, double x2, int pt, int n, int i, int crossing){
	int N = sweeps.xvlist[i].N;
	
	if (N+2>sweeps.xvlist[i]._N){
		sweeps.xvlist[i]._N 		+= 1024;
		sweeps.xvlist[i].xvalues	= (struct xvalue*)realloc(sweeps.xvlist[i].xvalues,sweeps.xvlist[i]._N*sizeof(struct xvalue));
		printf("xvlist size now %d\n",sweeps.xvlist[i]._N);
	}

	sweeps.xvlist[i].xvalues[N].x 	= x1;
	sweeps.xvlist[i].xvalues[N].pt 	= pt;
	sweeps.xvlist[i].xvalues[N].nx 	= n;
	sweeps.xvlist[i].xvalues[N].inout 	= 0;
	sweeps.xvlist[i].xvalues[N].crossing 	= crossing;
	sweeps.xvlist[i].xvalues[N+1].x 	= x2;
	sweeps.xvlist[i].xvalues[N+1].pt 	= pt;
	sweeps.xvlist[i].xvalues[N+1].nx 	= n;
	sweeps.xvlist[i].xvalues[N+1].inout	= 1;
	sweeps.xvlist[i].xvalues[N+1].crossing= crossing;

	sweeps.xvlist[i].N += 2;
}

void add_line_to_xvlist(double x1, double x2, int pt, int n, int crossing){
	int ix1 = (int)(floor( (x1/sweeps.widthx+0.5) *(double)sweeps.proc));// %sweeps.xvlists;
	int ix2 = (int)(floor( (x2/sweeps.widthx+0.5) *(double)sweeps.proc));// %sweeps.xvlists;
	if (ix2>=sweeps.proc){
		ix2 = sweeps.proc-1;
	}
	if (ix1<0){
		ix1 = 0;
	}

	if (ix1!=ix2){
		double b = -sweeps.widthx/2.+sweeps.widthx/(double)sweeps.proc*(double)ix2; 
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

	if (xmin<-sweeps.widthx/2.){
		add_line_to_xvlist(xmin+sweeps.widthx,sweeps.widthx/2.,pt,1,1);
		add_line_to_xvlist(-sweeps.widthx/2.,xmax,pt,0,1);
		return;
	}
	if (xmax>sweeps.widthx/2.){
		add_line_to_xvlist(-sweeps.widthx/2.,xmax-sweeps.widthx,pt,-1,1);
		add_line_to_xvlist(xmin,sweeps.widthx/2.,pt,0,1);
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
	
		double oldx = particles[i].x-dt*particles[i].vx;	
		add_to_xvlist(oldx,particles[i].x,i);
	}
//#pragma omp parallel for
	for (int proc=0;proc<sweeps.proc;proc++){
		struct xvaluelist xvlist = sweeps.xvlist[proc];
		qsort (xvlist.xvalues, xvlist.N, sizeof(struct xvalue), compare_xvalue);

		struct xvalue** sweeps 	= malloc(sizeof(struct xvalue*)*1024); // Active list. Make the magic number go away.
		int sweeps_N			= 0;

		for (int i=0;i<xvlist.N;i++){
			struct xvalue* const xv = &(xvlist.xvalues[i]);
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
					int gbnz = 0;
					for (int gbny = -1; gbny<=1; gbny++){
						struct ghostbox gb = get_ghostbox(gbnx,gbny,gbnz);
						detect_collision_of_pair(p1,p2,gb,proc,sweeps[k]->crossing||xv->crossing);
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

double min(double a, double b){ return (a>b)?b:a;}

double max(double a, double b){ return (b>a)?b:a;}

double sgn(const double a){
	return (a>=0 ? 1. : -1);
}

void detect_collision_of_pair(const int pt1, const int pt2, struct ghostbox const gb, int proc, int crossing){
	struct particle* p1 = &(particles[pt1]);
	struct particle* p2 = &(particles[pt2]);
	const double y 	= p1->y	+ gb.shifty	- p2->y;
	const double vy = p1->vy	+ gb.shiftvy 	- p2->ovy;

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
			struct collisionlist* const clist = &(sweeps.clist[proc]);
			if (clist->N>=clist->_N){
				if (clist->_N==0){
					clist->_N 		= 1024;
					clist->collisions	= (struct collision*)malloc(clist->_N*sizeof(struct collision));
				}else{
					clist->_N	 	+= 1024;
					clist->collisions	= (struct collision*)realloc(clist->collisions,clist->_N*sizeof(struct collision));
				}
				printf("Collisions array size now: %d\n", clist->_N);
			}
			struct collision* const c = &(clist->collisions[clist->N]);
			c->p1		= pt1;
			c->p2		= pt2;
			c->gb	 	= gb;
			c->time 	= timesave;
			c->crossing 	= crossing;
			clist->N++;
		}
	}
}

int compare_time(const void * a, const void * b){
	struct collision* ca = (struct collision*)a;
	struct collision* cb = (struct collision*)b;
	const double diff = ca->time - cb->time;
	if (diff > 0) return 1;
	if (diff < 0) return -1;
	return 0;
}
void collisions_resolve(){
#ifdef _OPENMP
	omp_lock_t boundarylock;
	omp_init_lock(&boundarylock);
#endif //_OPENMP

//#pragma omp parallel for
	for (int proc=0;proc<sweeps.proc;proc++){
		struct collision* c = sweeps.clist[proc].collisions;
		int N = sweeps.clist[proc].N;
	
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
			collisions_resolve_single(c1);
#ifdef _OPENMP
			if (c1.crossing){
				omp_unset_lock(&boundarylock);
			}
#endif //_OPENMP
		}
		sweeps.clist[proc].N = 0;
		sweeps.xvlist[proc].N = 0;
	}
#ifdef _OPENMP
	omp_destroy_lock(&boundarylock);
#endif //_OPENMP
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


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

