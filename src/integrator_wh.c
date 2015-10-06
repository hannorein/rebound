/**
 * @file 	integrator_wh.c
 * @brief 	Wisdom-Holman integrator.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	This file implements the Wisdom-Holman integration scheme.  
 * This scheme is second order accurate, symplectic and well suited for 
 * systems where there is one dominant mass and all particles are nearly on 
 * Keplerian orbits. Note that the scheme is formally only first order 
 * accurate when velocity dependent forces are present.
 * reb_particles should be sorted by increasing semi-major axis because the 
 * integrator transforms the positions and velocities to Jacobi coordinates
 * during the timestep and assumes that the particles are sorted.
 * The code is based on SWIFT. 
 * References: AJ, vol. 102, Oct. 1991, p. 1528-1538.
 *
 * The central mass is fixed at the origin (x,y,z)=(0,0,0). The indirect term
 * can be included by setting WH_INDIRECT_TERM equal to 1. When self-gravity 
 * is not included, the code can be run with a single sub-step insted of two by
 * setting WH_SELF_GRAVITY_ENABLED equal to 1. 
 * 
 * @section 	LICENSE
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
#include "rebound.h"
#include "particle.h"
#include "gravity.h"
#include "integrator.h"
#include "integrator_wh.h"

static void reb_drift_wh(struct reb_particle* const particles, const double G, double _dt, const int N, const int N_active);
static void reb_drift_dan(struct reb_particle* pv, double mu, double dt, int* iflag);
static void reb_drift_kepu(double dt, double r0, double mu, double alpha, double u, double* fp, double* c1, double* c2, double* c3, int* iflag);
static void reb_drift_kepu_guess(double dt0, double r0, double mu, double alpha, double u, double* s);
static void reb_drift_kepu_p3solve(double dt0, double r0, double mu, double alpha, double u, double* s, int* iflag);
static void reb_drift_kepu_new(double* s, double dt0, double r0, double mu, double alpha, double u, double* fp, double* c1, double* c2, double* c3, int* iflag);
static void reb_drift_kepu_lag(double* s, double dt0, double r0, double mu, double alpha, double u, double* fp, double* c1, double* c2, double* c3, int* iflag);
static void reb_drift_kepu_stumpff(double x, double* c0, double* c1, double* c2, double* c3);
static void reb_drift_kepu_fchk(double dt0, double r0, double mu, double alpha, double u, double s, double* f);
static void reb_drift_kepmd(double dm, double es, double ec, double* x, double* s, double* c);
static void reb_integrator_wh_aj(struct reb_particle* const particles, const double G, const int N, const int N_active);
static void reb_integrator_wh_ah(struct reb_particle* const particles, const double G, const int N, const int N_active);
static void reb_integrator_wh_to_jacobi(struct reb_particle* const particles, const double* const eta, const int N, const int N_active);
static void reb_integrator_wh_from_jacobi(struct reb_particle* const particles, const double* const eta, const int N, const int N_active);

#define DANBYB 1.e-13	///< Close to smallest relative floating point number

void reb_integrator_wh_part1(struct reb_simulation* r){
	const int N = r->N;
	const int N_active = r->N_active;
	struct reb_particle* const particles = r->particles;
	int _N_active = (N_active==-1)?N:N_active;
	if (_N_active!=r->ri_wh.allocatedN){
		r->ri_wh.eta = realloc(r->ri_wh.eta,sizeof(double)*_N_active);
		r->ri_wh.allocatedN = _N_active;
	}
	// DRIFT
	double* const eta = r->ri_wh.eta;
	eta[0] = particles[0].m;
	for(int i=1;i<_N_active;i++){
	  eta[i] = eta[i-1] + particles[i].m;
	}
	reb_integrator_wh_to_jacobi(particles, eta, N, N_active);
	reb_drift_wh(particles, r->G, r->dt/2., N, N_active);
	reb_integrator_wh_from_jacobi(particles, eta, N, N_active);
	r->t+=r->dt/2.;
}

void reb_integrator_wh_part2(struct reb_simulation* r){
	const int N = r->N;
	const int N_active = r->N_active;
	struct reb_particle* const particles = r->particles;
	const double* const eta = r->ri_wh.eta;
	// KICK
	// Calculate terms in Heliocentric coordinates
	reb_integrator_wh_ah(particles, r->G, N, N_active);
	reb_integrator_wh_to_jacobi(particles, eta, N, N_active);
	// Calculate terms in Jacobi coordinates
	reb_integrator_wh_aj(particles, r->G, N, N_active);
	reb_integrator_wh_from_jacobi(particles, eta, N, N_active);
	for (int i=1;i<N;i++){
		particles[i].vx += r->dt*particles[i].ax;
		particles[i].vy += r->dt*particles[i].ay;
		particles[i].vz += r->dt*particles[i].az;
	}
	// DRIFT
	reb_integrator_wh_to_jacobi(particles, eta, N, N_active);
	reb_drift_wh(particles, r->G, r->dt/2., N, N_active);
	reb_integrator_wh_from_jacobi(particles, eta, N, N_active);
	r->t+=r->dt/2.;
	r->dt_last_done = r->dt;
}

void reb_integrator_wh_synchronize(struct reb_simulation* r){
	// Do nothing.
}
void reb_integrator_wh_reset(struct reb_simulation* r){
	free(r->ri_wh.eta);
	r->ri_wh.allocatedN = 0;
}


// Only static routines below

static int wh_check_normal(struct reb_particle* p){
	if (isnan(p->vx) || isnan(p->vy)) return 1;
	if (isinf(p->vx) || isinf(p->vy)) return 2;
	return 0;
}

static void reb_integrator_wh_to_jacobi(struct reb_particle* const particles, const double* const eta, const int N, const int N_active){
	int _N_active = (N_active==-1)?N:N_active;

	double sumx  = particles[1].m * particles[1].x;
	double sumy  = particles[1].m * particles[1].y;
	double sumz  = particles[1].m * particles[1].z;
	double sumvx = particles[1].m * particles[1].vx;
	double sumvy = particles[1].m * particles[1].vy;
	double sumvz = particles[1].m * particles[1].vz;

	double capx  = sumx/eta[1];
	double capy  = sumy/eta[1];
	double capz  = sumz/eta[1];
	double capvx = sumvx/eta[1];
	double capvy = sumvy/eta[1];
	double capvz = sumvz/eta[1];

	for (int i=2;i<_N_active;i++){
		struct reb_particle p = particles[i];
		sumx  += p.m*p.x;
		sumy  += p.m*p.y;
		sumz  += p.m*p.z;
		sumvx += p.m*p.vx;
		sumvy += p.m*p.vy;
		sumvz += p.m*p.vz;
		
		particles[i].x  -= capx;
		particles[i].y  -= capy;
		particles[i].z  -= capz;
		particles[i].vx -= capvx;
		particles[i].vy -= capvy;
		particles[i].vz -= capvz;

		capx  = sumx /eta[i];
		capy  = sumy /eta[i];
		capz  = sumz /eta[i];
		capvx = sumvx/eta[i];
		capvy = sumvy/eta[i];
		capvz = sumvz/eta[i];
	}
}

static void reb_integrator_wh_from_jacobi(struct reb_particle* const particles, const double* const eta, const int N, const int N_active){
	int _N_active = (N_active==-1)?N:N_active;
	double sumx  = particles[1].m*particles[1].x /eta[1];
	double sumy  = particles[1].m*particles[1].y /eta[1];
	double sumz  = particles[1].m*particles[1].z /eta[1];
	double sumvx = particles[1].m*particles[1].vx/eta[1];
	double sumvy = particles[1].m*particles[1].vy/eta[1];
	double sumvz = particles[1].m*particles[1].vz/eta[1];

	for(int i=2;i<_N_active;i++){
		struct reb_particle p = particles[i];

		particles[i].x  += sumx;
		particles[i].y  += sumy;
		particles[i].z  += sumz;
		particles[i].vx += sumvx;
		particles[i].vy += sumvy;
		particles[i].vz += sumvz;

		sumx  += p.m*p.x  / eta[i];
		sumy  += p.m*p.y  / eta[i];
		sumz  += p.m*p.z  / eta[i];
		sumvx += p.m*p.vx / eta[i];
		sumvy += p.m*p.vy / eta[i];
		sumvz += p.m*p.vz / eta[i];
	}
}

// Assumes positions in heliocentric coordinates
static void reb_integrator_wh_ah(struct reb_particle* const particles, const double G, const int N, const int N_active){
	int _N_active = (N_active==-1)?N:N_active;
	// Massive particles
	double mass0 = particles[0].m;
	double axh0 = 0.0;
	double ayh0 = 0.0;
	double azh0 = 0.0;
	for(int i=2;i<_N_active;i++){
		struct reb_particle p = particles[i];
		double r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
		double ir3h = 1./(r*r*r);
		double fac = G*p.m*ir3h;
		axh0 -= fac*p.x;
		ayh0 -= fac*p.y;
		azh0 -= fac*p.z;
	}
	for(int i=1;i<_N_active;i++){
		struct reb_particle* p = &(particles[i]);
		double r = sqrt(p->x*p->x + p->y*p->y + p->z*p->z);
		double ir3h = 1./(r*r*r);
		if (i==1) ir3h=0.;

		p->ax += axh0 - G*mass0*p->x*ir3h;
		p->ay += ayh0 - G*mass0*p->y*ir3h;
		p->az += azh0 - G*mass0*p->z*ir3h;
	}
	// Test particles
	struct reb_particle p = particles[1];
	double r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
	double ir3h = 1./(r*r*r);
	double fac = G*p.m*ir3h;
	axh0 -= fac*p.x;
	ayh0 -= fac*p.y;
	azh0 -= fac*p.z;

#pragma omp parallel for schedule(guided)
	for(int i=_N_active;i<N;i++){
		struct reb_particle* pt = &(particles[i]);

		pt->ax += axh0;
		pt->ay += ayh0;
		pt->az += azh0;
	}
}

// Assumes position in Jacobi coordinates
static void reb_integrator_wh_aj(struct reb_particle* const particles, const double G, const int N, const int N_active){
	int _N_active = (N_active==-1)?N:N_active;
	// Massive particles (No need to calculate this for test particles)
	double mass0 = particles[0].m;
	double etaj = mass0;
	double axh2 = 0;
	double ayh2 = 0;
	double azh2 = 0;
	for(int i=2;i<_N_active;i++){
		struct reb_particle* p = &(particles[i]);
		double rj = sqrt(p->x*p->x + p->y*p->y + p->z*p->z);
		double ir3j = 1./(rj*rj*rj);
		
		etaj = etaj + particles[i-1].m;
		double fac = G*p->m*mass0*ir3j/etaj;
		axh2 += fac*p->x;
		ayh2 += fac*p->y;
		azh2 += fac*p->z;

		p->ax += G*mass0*p->x*ir3j + axh2;
		p->ay += G*mass0*p->y*ir3j + ayh2;
		p->az += G*mass0*p->z*ir3j + azh2;
	}
}

/**
 * @brief This function integrates the Keplerian motion of all particles.
 */
static void reb_drift_wh(struct reb_particle* const particles, const double G, double _dt, const int N, const int N_active){
	int _N_active = (N_active==-1)?N:N_active;
	double mass0 = particles[0].m;
	double etajm1 = mass0;
	// Massive particles
	for (int i=1;i<_N_active;i++){
		struct reb_particle* p = &(particles[i]);
		double etaj = etajm1 + p->m;
		double mu = G*mass0*etaj/etajm1;
		if (wh_check_normal(p)!=0){
			etaj = etajm1;
			continue;
		}
		int iflag = 0;
		reb_drift_dan(p,mu,_dt,&iflag);
		if (iflag != 0){ // Try again with 10 times smaller timestep.
			for (int j=0;j<10;j++){
				reb_drift_dan(p,mu,_dt/10.,&iflag);
				if (iflag != 0) break;
			}
		}
		etajm1 = etaj;  // Fixed by Subo
	}
	// Testparticles
#pragma omp parallel for schedule(guided)
	for (int i=_N_active;i<N;i++){
		struct reb_particle* p = &(particles[i]);
		if (wh_check_normal(p)!=0) continue;
		int iflag = 0;
		reb_drift_dan(p,G*mass0,_dt,&iflag);
		if (iflag != 0){ // Try again with 10 times smaller timestep.
			for (int j=0;j<10;j++){
				reb_drift_dan(p,G*mass0,_dt/10.,&iflag);
				if (iflag != 0) break;
			}
		}
	}
}

static void reb_drift_dan(struct reb_particle* pv, double mu, double dt0, int* iflag){
	double dt1 = dt0;
	double x0 = pv->x;
	double y0 = pv->y;
	double z0 = pv->z;
	double vx0 = pv->vx;
	double vy0 = pv->vy;
	double vz0 = pv->vz;
	
	double r0 = sqrt(x0*x0 + y0*y0 + z0*z0);
	double v0s = vx0*vx0 + vy0*vy0 + vz0*vz0;
	double u = x0*vx0 + y0*vy0 + z0*vz0;
	double alpha = 2.0*mu/r0 - v0s;

	if (alpha > 0.0){
		double a = mu/alpha;
		double asq = a*a;
		double en = sqrt(mu/(a*asq));
		double ec = 1.0 - r0/a;
		double es = u/(en*asq);
		double esq = ec*ec + es*es;
		double dm = dt1*en - floor(dt1*en/(2.0*M_PI))*2.0*M_PI; // TODO: Check that floor = int in fortran
		dt1 = dm/en;
		if ((esq*dm*dm < 0.0016) && !(dm*dm > 0.16 || esq > 0.36) ){
			double s, c, xkep;
			reb_drift_kepmd(dm,es,ec,&xkep,&s,&c);
			double fchk = (xkep - ec*s + es*(1.-c) - dm);
			if (fchk*fchk > DANBYB){
				*iflag =1;
				return;
			}
			
			double fp = 1. - ec*c + es*s;
			double f = (a/r0) * (c-1.) + 1.;
			double g = dt1 + (s-xkep)/en;
			double fdot = - (a/(r0*fp))*en*s;
			double gdot = (c-1.)/fp + 1.;
	
			pv->x = x0*f + vx0*g;
			pv->y = y0*f + vy0*g;
			pv->z = z0*f + vz0*g;
			
			pv->vx = x0*fdot + vx0*gdot;
			pv->vy = y0*fdot + vy0*gdot;
			pv->vz = z0*fdot + vz0*gdot;

			*iflag =0;
			return;
		}
	}
	double c1;
	double c2;
	double c3;
	double fp;
	reb_drift_kepu(dt1,r0,mu,alpha,u,&fp,&c1,&c2,&c3,iflag);
	if (*iflag==0){
		double f = 1.0 - (mu/r0)*c2;
		double g = dt1 - mu*c3;
		double fdot = -(mu/(fp*r0))*c1;
		double gdot = 1.0 - (mu/fp)*c2;

		pv->x = x0*f + vx0*g;
		pv->y = y0*f + vy0*g;
		pv->z = z0*f + vz0*g;
		
		pv->vx = x0*fdot + vx0*gdot;
		pv->vy = y0*fdot + vy0*gdot;
		pv->vz = z0*fdot + vz0*gdot;
	}
}


static void reb_drift_kepu(double dt0, double r0, double mu, double alpha, double u,
		double* fp, double* c1, double* c2, double* c3, int* iflag){
	// iflag == 0 if converged
	double s, st;
	reb_drift_kepu_guess(dt0, r0, mu, alpha, u, &s);
	st = s;
	reb_drift_kepu_new(&s,dt0,r0,mu,alpha,u,fp,c1,c2,c3,iflag);
	if (*iflag!=0){
		// fall back
		double fo,fn;
		reb_drift_kepu_fchk(dt0,r0,mu,alpha,u,st,&fo);
		reb_drift_kepu_fchk(dt0,r0,mu,alpha,u,s, &fn);
		if (fabs(fo)<fabs(fn)){
			s = st;
		}
		reb_drift_kepu_lag(&s,dt0,r0,mu,alpha,u,fp,c1,c2,c3,iflag);
	}
}

static void reb_drift_kepu_guess(double dt0, double r0, double mu, double alpha, double u, double* s){
	if (alpha > 0.){
		// elliptic motion
		if ((dt0/r0) <= 0.4){
			*s = dt0/r0 - (dt0*dt0*u)/(2.*r0*r0*r0);
			return;
		}else{
			double a = mu/alpha;
			double en = sqrt(mu/(a*a*a));
			double ec = 1. - r0/a;
			double es = u/(en*a*a);
			double e = sqrt(ec*ec + es*es);
			double y = en*dt0 - es;
			double sy = sin(y);
			double cy = cos(y);
			double sigma = ( (es*cy + ec*sy)>=0 ? 1. : -1. ); 
			double x = y + sigma*0.85*e;
			*s = x/sqrt(alpha);
		}
	}else{
		// hyperbolic
		int iflag=0;
		reb_drift_kepu_p3solve(dt0,r0,mu,alpha,u,s,&iflag);
		if (iflag!=0){
			*s = dt0/r0;
		}
	}
}

static void reb_drift_kepu_p3solve(double dt0, double r0, double mu, double alpha, double u, double* s, int* iflag){
	double denom = (mu - alpha*r0)/6.;
	double a2 = 0.5*u/denom;
	double a1 = r0/denom;
	double a0 = -dt0/denom;
	double q = (a1 - a2*a2/3.)/3.;
	double r = (a1*a2 - 3.*a0)/6. - a2*a2*a2/27.;
	double sq2 = q*q*q + r*r;

	if (sq2 >= 0.) {
		double sq = sqrt(sq2);
		double p1, p2;
		if ((r+sq) <= 0.){
			p1 = -pow(-(r+sq),1./3.);
		}else{
			p1 = pow(r+sq,1./3.);
		}
		if ((r-sq) <= 0.){
			p2 = -pow(-(r-sq),1./3.);
		}else{
			p2 = pow(r-sq,1./3.);
		}
		*iflag = 0;
		*s = p1 + p2 - a2/3.;
	}else{
		*iflag = 1;
		*s = 0;
	}
}

static void reb_drift_kepu_new(double* s, double dt0, double r0, double mu, double alpha, double u, double* fp, double* c1, double* c2, double* c3, int* iflag){
	int nc;
	for (nc=0;nc<6;nc++){
		double x = (*s)*(*s)*alpha;
		double c0;
		reb_drift_kepu_stumpff(x,&c0,c1,c2,c3);
		*c1 = (*c1) * (*s);
		*c2 = (*c2) * (*s) * (*s);
		*c3 = (*c3) * (*s) * (*s) * (*s);
		double f = r0*(*c1) + u*(*c2) + mu*(*c3) - dt0;
		(*fp) = r0*c0 + u*(*c1) + mu*(*c2);
		double fpp = (-r0*alpha + mu)*(*c1) + u*c0;
		double fppp = (-r0*alpha + mu)*c0 - u*alpha*(*c1);
		double ds = -f/(*fp);
		ds = -f/((*fp) +ds*fpp/2.);
		ds = -f/((*fp) +ds*fpp/2.+ds*ds*fppp/6.);
		*s = (*s) + ds;
		double fdt = f/dt0;
		if (fdt*fdt < DANBYB*DANBYB){
			*iflag = 0;
			return;
		}
	}
	// Not converged
	*iflag = 1;
}

static void reb_drift_kepu_stumpff(double x, double* c0, double* c1, double* c2, double* c3){
	int n = 0;
	double xm = 0.1;
	while(fabs(x)>= xm){
		n++;
		x = x / 4.;
	}
	
	*c2 = (1.-x*(1.-x*(1.-x*(1.-x*(1.-x*(1.-x/182.)/132.)/90.)/56.)/30.)/12.)/2.;
	*c3 = (1.-x*(1.-x*(1.-x*(1.-x*(1.-x*(1.-x/210.)/156.)/110.)/72.)/42.)/20.)/6.;
	*c1 = 1. - x*(*c3);
	*c0 = 1. - x*(*c2);
	for (int i=n;i>=1;i--){
		*c3 = ((*c2) + (*c0)*(*c3))/4.;
		*c2 = (*c1)*(*c1)/2.;
		*c1 = (*c0)*(*c1);
		*c0 = 2.*(*c0)*(*c0) - 1.;
		x = x * 4.; 
	}
}

static void reb_drift_kepu_fchk(double dt0, double r0, double mu, double alpha, double u, double s, double* f){
	double x = s*s*alpha;
	double c0, c1, c2, c3;
	reb_drift_kepu_stumpff(x,&c0,&c1,&c2,&c3);
	c1 = c1 *s;
	c2 = c2 *s*s;
	c3 = c3 *s*s*s;
	*f = r0*c1 + u*c2 + mu*c3 - dt0;
}

static void reb_drift_kepu_lag(double* s, double dt0, double r0, double mu, double alpha, double u, double* fp, double* c1, double* c2, double* c3, int* iflag){
	int ncmax = 400;
	double ln = 5.;
	
	int nc;
	for (nc=0;nc<=ncmax;nc++){
		double x = (*s)*(*s)*alpha;
		double c0;
		reb_drift_kepu_stumpff(x,&c0,c1,c2,c3);
		*c1 = (*c1)*(*s);
		*c2 = (*c2)*(*s)*(*s);
		*c3 = (*c3)*(*s)*(*s)*(*s);
		double f = r0*(*c1) + u*(*c2) + mu*(*c3) - dt0;
		(*fp) = r0*c0 + u*(*c1) + mu*(*c2);
		double fpp = (-40.*alpha + mu)*(*c1) + u*c0;
		double ds = -ln*f/((*fp) + ((*fp)>0.?1.:-1.)*sqrt(fabs((ln - 1.) * (ln - 1.) *(*fp)*(*fp) - (ln - 1.) * ln * f * fpp)));
		*s = (*s) + ds;

		double fdt = f/dt0;
		if (fdt*fdt < DANBYB*DANBYB){
			*iflag = 0;
			return;
		}
	}
	*iflag = 2;
}

static void reb_drift_kepmd(double dm, double es, double ec, double* x, double* s, double* c){
	const double A0 = 39916800.;
	const double A1 = 6652800.;
	const double A2 = 332640.;
	const double A3 = 7920.;
	const double A4 = 110.;

	double fac1 = 1./(1. -ec);
	double q = fac1 * dm;
	double fac2 = es*es*fac1 - ec/3.;
	*x = q*(1. - 0.5*fac1*q*(es-q*fac2));

	double y = (*x)*(*x);
	*s = (*x)*(A0-y*(A1-y*(A2-y*(A3-y*(A4-y)))))/A0;
	*c = sqrt(1. - (*s)*(*s));

	double f = (*x) - ec*(*s) + es*(1.-(*c)) -dm;
	double fp = 1. - ec*(*c) + es*(*s);
	double fpp = ec*(*s) + es*(*c);
	double fppp = ec*(*c) - es*(*s);
	double dx = -f/fp;
	dx = -f/(fp + 0.5*dx*fpp);
	dx = -f/(fp + 0.5*dx*fpp + 0.16666666666666666666*dx*dx*fppp);
	*x = (*x) + dx;

	y = (*x)*(*x);
	*s = (*x)*(A0-y*(A1-y*(A2-y*(A3-y*(A4-y)))))/A0;
	*c = sqrt(1. - (*s)*(*s));
}

