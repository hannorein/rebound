/**
 * @file 	integrator.c
 * @brief 	Mikkola integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the WHFast integration scheme.  
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Hanno Rein, Daniel Tamayo
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
#include <sys/time.h>
#include "particle.h"
#include "main.h"
#include "tools.h"
#include "gravity.h"
#include "boundaries.h"
#include "integrator.h"
#include "integrator_whfast.h"

#define MAX(a, b) ((a) < (b) ? (b) : (a))
#define MIN(a, b) ((a) > (b) ? (b) : (a))

// the defaults below are chosen to safeguard the user against spurious results, but
// will be slower and less accurate
unsigned int integrator_whfast_safe_mode = 1;
unsigned int integrator_whfast_corrector = 0;
unsigned int integrator_whfast_recalculate_jacobi_this_timestep	= 0;

static unsigned int integrator_whfast_is_synchronized = 1;
static unsigned int integrator_allocated_N = 0;
static struct particle* restrict p_j  = NULL;
static double* restrict eta = NULL;
static double* restrict etai = NULL;
static double Mtotali;
static double Mtotal;
static unsigned int integrator_timestep_warning = 0;
static unsigned int integrator_whfast_recalculate_jacobi_but_not_synchronized_warning = 0;

// Fast inverse factorial lookup table
static const double invfactorial[35] = {1., 1., 1./2., 1./6., 1./24., 1./120., 1./720., 1./5040., 1./40320., 1./362880., 1./3628800., 1./39916800., 1./479001600., 1./6227020800., 1./87178291200., 1./1307674368000., 1./20922789888000., 1./355687428096000., 1./6402373705728000., 1./121645100408832000., 1./2432902008176640000., 1./51090942171709440000., 1./1124000727777607680000., 1./25852016738884976640000., 1./620448401733239439360000., 1./15511210043330985984000000., 1./403291461126605635584000000., 1./10888869450418352160768000000., 1./304888344611713860501504000000., 1./8841761993739701954543616000000., 1./265252859812191058636308480000000., 1./8222838654177922817725562880000000., 1./263130836933693530167218012160000000., 1./8683317618811886495518194401280000000., 1./295232799039604140847618609643520000000.};

//static double ipow(double base, unsigned int exp) {
//	double result = 1;
//	while (exp) {
//		if (exp & 1)
//		    result *= base;
//		exp >>= 1;
//		base *= base;
//	}
//	return result;
//}

static inline double fastabs(double x){
	    return (x > 0.) ? x : -x;
}

static void stumpff_cs(double *restrict cs, double z) {
	unsigned int n = 0;
	while(z>0.1){
		z = z/4.;
		n++;
	}
	double zm = -z;
	cs[4] = invfactorial[4] - z*invfactorial[6]; 	// always calculate first two terms
	cs[5] = invfactorial[5] - z*invfactorial[7]; 	// always calculate first two terms
	double old_c_4;
	double _pow = zm;
	unsigned int k=8;
	do{
		old_c_4 = cs[4];
		_pow *= zm;
		cs[4] += _pow*invfactorial[k];
		k+=1;
		cs[5] += _pow*invfactorial[k];
		k+=1;
	}while(cs[4]!=old_c_4 && k<35);			// Stop if new term smaller than machine precision (cs[5] converges faster than cs[4])
	cs[3] = 1./6.-z*cs[5];
	cs[2] = 0.5-z*cs[4];
	cs[1] = 1.-z*cs[3];
	for (;n>0;n--){	
		z = z*4.;
		cs[5] = (cs[5]+cs[4]+cs[3]*cs[2])*0.0625;
		cs[4] = (1.+cs[1])*cs[3]*0.125;
		cs[3] = 1./6.-z*cs[5];
		cs[2] = 0.5-z*cs[4];
		cs[1] = 1.-z*cs[3];
	}
	cs[0] = 1.-z*cs[2];
}

static void stumpff_cs3(double *restrict cs, double z) {
	unsigned int n = 0;
	while(z>0.1){
		z = z/4.;
		n++;
	}
	double zm = -z;
	cs[2] = invfactorial[2] - z*invfactorial[4]; 	// always calculate first two terms
	cs[3] = invfactorial[3] - z*invfactorial[5]; 	// always calculate first two terms
	double old_c_2;
	double _pow = zm;
	unsigned int k=6;
	do{
		old_c_2 = cs[2];
		_pow *= zm;
		cs[2] += _pow*invfactorial[k];
		k+=1;
		cs[3] += _pow*invfactorial[k];
		k+=1;
	}while(cs[2]!=old_c_2 && k<34);			// Stop if new term smaller than machine precision (cs[3] converges faster than cs[2])
	cs[1] = 1.-z*cs[3];
	cs[0] = 1.-z*cs[2];
	for (;n>0;n--){	
		cs[3] = (cs[2]+cs[0]*cs[3])*0.25;
		cs[2] = cs[1]*cs[1]*0.5;
		cs[1] = cs[0]*cs[1];
		cs[0] = 2.*cs[0]*cs[0]-1.;
	}
}

static void stiefel_Gs(double *restrict Gs, double beta, double X) {
	double X2 = X*X;
	stumpff_cs(Gs, beta*X2);
	Gs[1] *= X; 
	Gs[2] *= X2; 
	double _pow = X2*X;
	Gs[3] *= _pow; 
	_pow *= X;
	Gs[4] *= _pow; 
	_pow *= X;
	Gs[5] *= _pow; 
	return;
}

static void stiefel_Gs3(double *restrict Gs, double beta, double X) {
	double X2 = X*X;
	stumpff_cs3(Gs, beta*X2);
	Gs[1] *= X; 
	Gs[2] *= X2; 
	Gs[3] *= X2*X;
	return;
}

static inline double _M(unsigned int i){
  	return G*(eta[i]); // Hanno 1
	//return G*(eta[i-1]); // Hanno2 
	//return G*(eta[i-1]*particles[i].m*eta[i-1]/eta[i]/(eta[i-1]+particles[i].m*eta[i-1]/eta[i])); // reduced mass jacobi
	//return G*(eta[i-1]*particles[i].m/(eta[i-1]+particles[i].m)); // reduced mass
	//return G*(eta[i]/eta[i-1]*particles[0].m);   // SSD
}

/****************************** 
 * Keplerian motion           */
int iter;

static void kepler_step(unsigned int i,double _dt){
	const double M = _M(i);
	const struct particle p1 = p_j[i];

	const double r0 = sqrt(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z);
	const double r0i = 1./r0;
	const double v2 =  p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz;
	const double beta = 2.*M*r0i - v2;
	const double eta0 = p1.x*p1.vx + p1.y*p1.vy + p1.z*p1.vz;
	const double zeta0 = M - beta*r0;
double X;
	double Gs[6]; 
		
	if (beta>0.){
		// Elliptic orbit
		//X = _dt*invperiod*X_per_period; // first order guess 
		double dtr0i = _dt*r0i;
		//X = dtr0i; // first order guess
		X = dtr0i * (1. - dtr0i*eta0*0.5*r0i); // second order guess
		//X = dtr0i *(1.- 0.5*dtr0i*r0i*(eta0-dtr0i*(eta0*eta0*r0i-1./3.*zeta0))); // third order guess
		//X = _dt*beta/M + eta0/M*(0.85*sqrt(1.+zeta0*zeta0/beta/eta0/eta0) - 1.);  // Dan's version 

	}else{
		// Hyperbolic orbit
		X = 0.; // Initial guess 
	}

	unsigned int converged = 0;
	double oldX = X; 

	stiefel_Gs3(Gs, beta, X);
	const double eta0Gs1zeta0Gs2 = eta0*Gs[1] + zeta0*Gs[2];
	double ri = 1./(r0 + eta0Gs1zeta0Gs2);
	X  = ri*(X*eta0Gs1zeta0Gs2-eta0*Gs[2]-zeta0*Gs[3]+_dt);
	double X_per_period = 2.*M_PI/sqrt(beta);

	if(fastabs(X-oldX) > 0.01*X_per_period){
		// Linear guess
		X = beta*_dt/M;
#define WHFAST_NMAX_QUART 64
		static double prevX[WHFAST_NMAX_QUART+1];
		for(int n_lag=1; n_lag < WHFAST_NMAX_QUART; n_lag++){
			stiefel_Gs3(Gs, beta, X);
			const double f = r0*X + eta0*Gs[2] + zeta0*Gs[3] - _dt;
			const double fp = r0 + eta0*Gs[1] + zeta0*Gs[2];
			const double fpp = eta0*Gs[0] + zeta0*Gs[1];
			const double denom = fp + sqrt(fabs(16.*fp*fp - 20.*f*fpp));
			X = (X*denom - 5.*f)/denom;
			iter++;	// DEBUG
			
			for(int i=1;i<n_lag;i++){
				if(X==prevX[i]){
					// Converged. Exit.
					n_lag = WHFAST_NMAX_QUART;
					converged = 1;
					break;
				}
			}
			prevX[n_lag] = X;
			
		}
		const double eta0Gs1zeta0Gs2 = eta0*Gs[1] + zeta0*Gs[2];
		ri = 1./(r0 + eta0Gs1zeta0Gs2);
	}else{
		double oldX2 = NAN; // NAN might be a GNU extension, any value other than X works.
#define WHFAST_NMAX_NEWT 32
		for (int n_hg=1;n_hg<WHFAST_NMAX_NEWT;n_hg++){
			oldX2 = oldX;
			oldX = X;
			stiefel_Gs3(Gs, beta, X);
			const double eta0Gs1zeta0Gs2 = eta0*Gs[1] + zeta0*Gs[2];
			ri = 1./(r0 + eta0Gs1zeta0Gs2);
			X  = ri*(X*eta0Gs1zeta0Gs2-eta0*Gs[2]-zeta0*Gs[3]+_dt);
			iter++;	// DEBUG
			
			if (X==oldX||X==oldX2){
				// Converged. Exit.
				converged = 1;
				break; 
			}
		}
	}

	if (converged == 0){ // Fallback to bisection 
		double X_min, X_max;
		if (beta>0.){
			//Elliptic
			double sqrt_beta = sqrt(beta);
			double invperiod = sqrt_beta*beta/(2.*M_PI*M);
			double X_per_period = 2.*M_PI/sqrt_beta;
			if (fabs(_dt)*invperiod>1. && integrator_timestep_warning == 0){
				integrator_timestep_warning++;
				fprintf(stderr,"\n\033[1mWarning!\033[0m Timestep is larger than at least one orbital period.\n");
			}
			X_min = X_per_period * floor(_dt*invperiod);
			X_max = X_min + X_per_period;
		}else{
			//Hyperbolic
			double h2 = r0*r0*v2-eta0*eta0;
			double q = h2/M/(1.+sqrt(1.-h2*beta/(M*M)));
			double vq = sqrt(h2)/q;
			X_min = 1./(vq+r0/_dt);
			X_max = _dt/q;
		}
		X = (X_max + X_min)/2.;
		do{
			iter++;	// DEBUG
			stiefel_Gs3(Gs, beta, X);
			double s   = r0*X + eta0*Gs[2] + zeta0*Gs[3]-_dt;
			if (s>=0.){
				X_max = X;
			}else{
				X_min = X;
			}
			X = (X_max + X_min)/2.;
		}while (fastabs((X_max-X_min)/X_max)>1e-15);
		const double eta0Gs1zeta0Gs2 = eta0*Gs[1] + zeta0*Gs[2];
		ri = 1./(r0 + eta0Gs1zeta0Gs2);
	}

	// Note: These are not the traditional f and g functions.
	double f = -M*Gs[2]*r0i;
	double g = _dt - M*Gs[3];
	double fd = -M*Gs[1]*r0i*ri; 
	double gd = -M*Gs[2]*ri; 
        
	p_j[i].x += f*p1.x + g*p1.vx;
	p_j[i].y += f*p1.y + g*p1.vy;
	p_j[i].z += f*p1.z + g*p1.vz;
        
	p_j[i].vx += fd*p1.x + gd*p1.vx;
	p_j[i].vy += fd*p1.y + gd*p1.vy;
	p_j[i].vz += fd*p1.z + gd*p1.vz;

	//Variations
	if (N_megno){
		stiefel_Gs(Gs, beta, X);	// Recalculate (to get Gs[4] and Gs[5])
		struct particle dp1 = p_j[i+N_megno];
		double dr0 = (dp1.x*p1.x + dp1.y*p1.y + dp1.z*p1.z)*r0i;
		double dbeta = -2.*M*dr0*r0i*r0i - 2.* (dp1.vx*p1.vx + dp1.vy*p1.vy + dp1.vz*p1.vz);
		double deta0 = dp1.x*p1.vx + dp1.y*p1.vy + dp1.z*p1.vz
			     + p1.x*dp1.vx + p1.y*dp1.vy + p1.z*dp1.vz;
		double dzeta0 = -beta*dr0 - r0*dbeta;
		double G3beta = 0.5*(3.*Gs[5]-X*Gs[4]);
		double G2beta = 0.5*(2.*Gs[4]-X*Gs[3]);
		double G1beta = 0.5*(Gs[3]-X*Gs[2]);
		double tbeta = eta0*G2beta + zeta0*G3beta;
		double dX = -1.*ri*(X*dr0 + Gs[2]*deta0+Gs[3]*dzeta0+tbeta*dbeta);
		double dG1 = Gs[0]*dX + G1beta*dbeta; 
		double dG2 = Gs[1]*dX + G2beta*dbeta;
		double dG3 = Gs[2]*dX + G3beta*dbeta;
		double dr = dr0 + Gs[1]*deta0 + Gs[2]*dzeta0 + eta0*dG1 + zeta0*dG2;
		double df = M*Gs[2]*dr0*r0i*r0i - M*dG2*r0i;
		double dg = -M*dG3;
		double dfd = -M*dG1*r0i*ri + M*Gs[1]*(dr0*r0i+dr*ri)*r0i*ri;
		double dgd = -M*dG2*ri + M*Gs[2]*dr*ri*ri;
	
		p_j[i+N_megno].x += f*dp1.x + g*dp1.vx + df*p1.x + dg*p1.vx;
		p_j[i+N_megno].y += f*dp1.y + g*dp1.vy + df*p1.y + dg*p1.vy;
		p_j[i+N_megno].z += f*dp1.z + g*dp1.vz + df*p1.z + dg*p1.vz;

		p_j[i+N_megno].vx += fd*dp1.x + gd*dp1.vx + dfd*p1.x + dgd*p1.vx;
		p_j[i+N_megno].vy += fd*dp1.y + gd*dp1.vy + dfd*p1.y + dgd*p1.vy;
		p_j[i+N_megno].vz += fd*dp1.z + gd*dp1.vz + dfd*p1.z + dgd*p1.vz;
	}

}

/****************************** 
 * Coordinate transformations */
static void integrator_to_jacobi_posvel(void){
	double s_x = particles[0].m * particles[0].x;
	double s_y = particles[0].m * particles[0].y;
	double s_z = particles[0].m * particles[0].z;
	double s_vx = particles[0].m * particles[0].vx;
	double s_vy = particles[0].m * particles[0].vy;
	double s_vz = particles[0].m * particles[0].vz;
	for (unsigned int i=1;i<N-N_megno;i++){
		const double ei = etai[i-1];
		const struct particle pi = particles[i];
		const double pme = eta[i]*ei;
		p_j[i].x = pi.x - s_x*ei;
		p_j[i].y = pi.y - s_y*ei;
		p_j[i].z = pi.z - s_z*ei;
		p_j[i].vx = pi.vx - s_vx*ei;
		p_j[i].vy = pi.vy - s_vy*ei;
		p_j[i].vz = pi.vz - s_vz*ei;
		s_x  = s_x  * pme + pi.m*p_j[i].x ;
		s_y  = s_y  * pme + pi.m*p_j[i].y ;
		s_z  = s_z  * pme + pi.m*p_j[i].z ;
		s_vx = s_vx * pme + pi.m*p_j[i].vx;
		s_vy = s_vy * pme + pi.m*p_j[i].vy;
		s_vz = s_vz * pme + pi.m*p_j[i].vz;
	}
	p_j[0].x = s_x * Mtotali;
	p_j[0].y = s_y * Mtotali;
	p_j[0].z = s_z * Mtotali;
	p_j[0].vx = s_vx * Mtotali;
	p_j[0].vy = s_vy * Mtotali;
	p_j[0].vz = s_vz * Mtotali;
}

static void integrator_var_to_jacobi_posvel(void){
	double s_x = particles[N_megno].m * particles[N_megno].x;
	double s_y = particles[N_megno].m * particles[N_megno].y;
	double s_z = particles[N_megno].m * particles[N_megno].z;
	double s_vx = particles[N_megno].m * particles[N_megno].vx;
	double s_vy = particles[N_megno].m * particles[N_megno].vy;
	double s_vz = particles[N_megno].m * particles[N_megno].vz;
	for (unsigned int i=1+N_megno;i<N;i++){
		const double ei = etai[i-1-N_megno];
		const struct particle pi = particles[i];
		const double pme = eta[i-N_megno]*ei;
		p_j[i].x = pi.x - s_x*ei;
		p_j[i].y = pi.y - s_y*ei;
		p_j[i].z = pi.z - s_z*ei;
		p_j[i].vx = pi.vx - s_vx*ei;
		p_j[i].vy = pi.vy - s_vy*ei;
		p_j[i].vz = pi.vz - s_vz*ei;
		s_x  = s_x  * pme + pi.m*p_j[i].x ;
		s_y  = s_y  * pme + pi.m*p_j[i].y ;
		s_z  = s_z  * pme + pi.m*p_j[i].z ;
		s_vx = s_vx * pme + pi.m*p_j[i].vx;
		s_vy = s_vy * pme + pi.m*p_j[i].vy;
		s_vz = s_vz * pme + pi.m*p_j[i].vz;
	}
	p_j[N_megno].x = s_x * Mtotali;
	p_j[N_megno].y = s_y * Mtotali;
	p_j[N_megno].z = s_z * Mtotali;
	p_j[N_megno].vx = s_vx * Mtotali;
	p_j[N_megno].vy = s_vy * Mtotali;
	p_j[N_megno].vz = s_vz * Mtotali;
}


static void integrator_to_jacobi_acc(void){
	double s_ax = particles[0].m * particles[0].ax;
	double s_ay = particles[0].m * particles[0].ay;
	double s_az = particles[0].m * particles[0].az;
	for (unsigned int i=1;i<N-N_megno;i++){
		const double ei = etai[i-1];
		const struct particle pi = particles[i];
		const double pme = eta[i]*ei;
		p_j[i].ax = pi.ax - s_ax*ei;
		p_j[i].ay = pi.ay - s_ay*ei;
		p_j[i].az = pi.az - s_az*ei;
		s_ax = s_ax * pme + pi.m*p_j[i].ax;
		s_ay = s_ay * pme + pi.m*p_j[i].ay;
		s_az = s_az * pme + pi.m*p_j[i].az;
	}
	// p_j[0].a is not needed and thus not calculated 
}

static void integrator_var_to_jacobi_acc(void){
	double s_ax = particles[N_megno].m * particles[N_megno].ax;
	double s_ay = particles[N_megno].m * particles[N_megno].ay;
	double s_az = particles[N_megno].m * particles[N_megno].az;
	for (unsigned int i=1+N_megno;i<N;i++){
		const double ei = etai[i-1-N_megno];
		const struct particle pi = particles[i];
		const double pme = eta[i-N_megno]*ei;
		p_j[i].ax = pi.ax - s_ax*ei;
		p_j[i].ay = pi.ay - s_ay*ei;
		p_j[i].az = pi.az - s_az*ei;
		s_ax = s_ax * pme + pi.m*p_j[i].ax;
		s_ay = s_ay * pme + pi.m*p_j[i].ay;
		s_az = s_az * pme + pi.m*p_j[i].az;
	}
	// p_j[N_megno].a is not needed and thus not calculated 
}

static void integrator_to_inertial_posvel(void){
	double s_x  = p_j[0].x  * Mtotal; 
	double s_y  = p_j[0].y  * Mtotal; 
	double s_z  = p_j[0].z  * Mtotal; 
	double s_vx = p_j[0].vx * Mtotal; 
	double s_vy = p_j[0].vy * Mtotal; 
	double s_vz = p_j[0].vz * Mtotal; 
	for (unsigned int i=N-N_megno-1;i>0;i--){
		const struct particle pji = p_j[i];
		const double ei = etai[i];
		s_x  = (s_x  - pji.m * pji.x ) * ei;
		s_y  = (s_y  - pji.m * pji.y ) * ei;
		s_z  = (s_z  - pji.m * pji.z ) * ei;
		s_vx = (s_vx - pji.m * pji.vx) * ei;
		s_vy = (s_vy - pji.m * pji.vy) * ei;
		s_vz = (s_vz - pji.m * pji.vz) * ei;
		particles[i].x  = pji.x  + s_x ;
		particles[i].y  = pji.y  + s_y ;
		particles[i].z  = pji.z  + s_z ;
		particles[i].vx = pji.vx + s_vx;
		particles[i].vy = pji.vy + s_vy;
		particles[i].vz = pji.vz + s_vz;
		s_x  *= eta[i-1];
		s_y  *= eta[i-1];
		s_z  *= eta[i-1];
		s_vx *= eta[i-1];
		s_vy *= eta[i-1];
		s_vz *= eta[i-1];
	}
	const double mi = etai[0];
	particles[0].x  = s_x  * mi;
	particles[0].y  = s_y  * mi;
	particles[0].z  = s_z  * mi;
	particles[0].vx = s_vx * mi;
	particles[0].vy = s_vy * mi;
	particles[0].vz = s_vz * mi;
}

static void integrator_var_to_inertial_posvel(void){
	double s_x  = p_j[N_megno].x  * Mtotal; 
	double s_y  = p_j[N_megno].y  * Mtotal; 
	double s_z  = p_j[N_megno].z  * Mtotal; 
	double s_vx = p_j[N_megno].vx * Mtotal; 
	double s_vy = p_j[N_megno].vy * Mtotal; 
	double s_vz = p_j[N_megno].vz * Mtotal; 
	for (unsigned int i=N-1;i>N_megno;i--){
		const struct particle pji = p_j[i];
		const double ei = etai[i-N_megno];
		s_x  = (s_x  - pji.m * pji.x ) * ei;
		s_y  = (s_y  - pji.m * pji.y ) * ei;
		s_z  = (s_z  - pji.m * pji.z ) * ei;
		s_vx = (s_vx - pji.m * pji.vx) * ei;
		s_vy = (s_vy - pji.m * pji.vy) * ei;
		s_vz = (s_vz - pji.m * pji.vz) * ei;
		particles[i].x  = pji.x  + s_x ;
		particles[i].y  = pji.y  + s_y ;
		particles[i].z  = pji.z  + s_z ;
		particles[i].vx = pji.vx + s_vx;
		particles[i].vy = pji.vy + s_vy;
		particles[i].vz = pji.vz + s_vz;
		s_x  *= eta[i-N_megno-1];
		s_y  *= eta[i-N_megno-1];
		s_z  *= eta[i-N_megno-1];
		s_vx *= eta[i-N_megno-1];
		s_vy *= eta[i-N_megno-1];
		s_vz *= eta[i-N_megno-1];
	}
	const double mi = etai[0];
	particles[N_megno].x  = s_x  * mi;
	particles[N_megno].y  = s_y  * mi;
	particles[N_megno].z  = s_z  * mi;
	particles[N_megno].vx = s_vx * mi;
	particles[N_megno].vy = s_vy * mi;
	particles[N_megno].vz = s_vz * mi;
}

static void integrator_to_inertial_pos(void){
	double s_x  = p_j[0].x  * Mtotal; 
	double s_y  = p_j[0].y  * Mtotal; 
	double s_z  = p_j[0].z  * Mtotal; 
	for (unsigned int i=N-N_megno-1;i>0;i--){
		const struct particle pji = p_j[i];
		const double ei = etai[i];
		s_x  = (s_x  - pji.m * pji.x ) * ei;
		s_y  = (s_y  - pji.m * pji.y ) * ei;
		s_z  = (s_z  - pji.m * pji.z ) * ei;
		particles[i].x  = pji.x  + s_x ;
		particles[i].y  = pji.y  + s_y ;
		particles[i].z  = pji.z  + s_z ;
		s_x  *= eta[i-1];
		s_y  *= eta[i-1];
		s_z  *= eta[i-1];
	}
	const double mi = etai[0];
	particles[0].x  = s_x  * mi;
	particles[0].y  = s_y  * mi;
	particles[0].z  = s_z  * mi;
}

static void integrator_var_to_inertial_pos(void){
	double s_x  = p_j[0].x  * Mtotal; 
	double s_y  = p_j[0].y  * Mtotal; 
	double s_z  = p_j[0].z  * Mtotal; 
	for (unsigned int i=N-1;i>N_megno;i--){
		const struct particle pji = p_j[i];
		const double ei = etai[i-N_megno];
		s_x  = (s_x  - pji.m * pji.x ) * ei;
		s_y  = (s_y  - pji.m * pji.y ) * ei;
		s_z  = (s_z  - pji.m * pji.z ) * ei;
		particles[i].x  = pji.x  + s_x ;
		particles[i].y  = pji.y  + s_y ;
		particles[i].z  = pji.z  + s_z ;
		s_x  *= eta[i-N_megno-1];
		s_y  *= eta[i-N_megno-1];
		s_z  *= eta[i-N_megno-1];
	}
	const double mi = etai[0];
	particles[N_megno].x  = s_x  * mi;
	particles[N_megno].y  = s_y  * mi;
	particles[N_megno].z  = s_z  * mi;
}

/***************************** 
 * Interaction Hamiltonian  */

static void integrator_interaction(double _dt){
	for (unsigned int i=1;i<N-N_megno;i++){
		// Eq 132
		const struct particle pji = p_j[i];
		static double rj2i;
		static double rj3iM;
		static double prefac1;
		p_j[i].vx += _dt * pji.ax;
		p_j[i].vy += _dt * pji.ay;
		p_j[i].vz += _dt * pji.az;
		if (i>1){
			rj2i = 1./(pji.x*pji.x + pji.y*pji.y + pji.z*pji.z + softening*softening);
			const double rji  = sqrt(rj2i);
			rj3iM = rji*rj2i*_M(i);
		 	prefac1 = _dt*rj3iM;
			p_j[i].vx += prefac1*pji.x;
			p_j[i].vy += prefac1*pji.y;
			p_j[i].vz += prefac1*pji.z;
		}
		if (N_megno){
			// Eq 132
			double rj5M = rj3iM*rj2i;
			double rdr = p_j[i+N_megno].x*pji.x + p_j[i+N_megno].y*pji.y + p_j[i+N_megno].z*pji.z;
			double prefac2 = -_dt*3.*rdr*rj5M;
			p_j[i+N_megno].vx += _dt * p_j[i+N_megno].ax;
			p_j[i+N_megno].vy += _dt * p_j[i+N_megno].ay;
			p_j[i+N_megno].vz += _dt * p_j[i+N_megno].az;
			if (i>1){
				p_j[i+N_megno].vx += prefac1*p_j[i+N_megno].x + prefac2*pji.x;
				p_j[i+N_megno].vy += prefac1*p_j[i+N_megno].y + prefac2*pji.y;
				p_j[i+N_megno].vz += prefac1*p_j[i+N_megno].z + prefac2*pji.z;
			}
		}
	}
}

/***************************** 
 * DKD Scheme                */

static void integrator_kepler_drift(double _dt){
	for (unsigned int i=1;i<N-N_megno;i++){
		kepler_step(i, _dt);
	}
	p_j[0].x += _dt*p_j[0].vx;
	p_j[0].y += _dt*p_j[0].vy;
	p_j[0].z += _dt*p_j[0].vz;
}

const static double a_1 = 4.183300132670377813e-01;
const static double a_2 = 2.*4.183300132670377813e-01;
const static double a_3 = 3.*4.183300132670377813e-01;
const static double a_4 = 4.*4.183300132670377813e-01;
const static double a_5 = 5.*4.183300132670377813e-01;
const static double b_31 = -0.5*4.980119205559973422e-02;
const static double b_51 = -1./6.*4.980119205559973422e-02;
const static double b_52 = 5./6.*4.980119205559973422e-02;
const static double b_71 = 12361./246960.*4.980119205559973422e-02;
const static double b_72 = -22651./61740.*4.980119205559973422e-02;
const static double b_73 = 53521./49392.*4.980119205559973422e-02;
const static double b_111 = 2798927./684573120.*4.980119205559973422e-02;
const static double b_112 = -329447./6985440.*4.980119205559973422e-02;
const static double b_113 = 895249./3622080.*4.980119205559973422e-02;
const static double b_114 = -14556229./19015920.*4.980119205559973422e-02;
const static double b_115 = 3394141./2328480.*4.980119205559973422e-02;


static void integrator_whfast_corrector_Z(double a, double b){
	integrator_kepler_drift(a);
	integrator_to_inertial_pos();
	if (N_megno){
		integrator_var_to_inertial_pos();
	}
	integrator_update_acceleration();
	integrator_to_jacobi_acc();
	if (N_megno){
		integrator_var_to_jacobi_acc();
	}
	integrator_interaction(-b);
	integrator_kepler_drift(-2.*a);
	integrator_to_inertial_pos();
	if (N_megno){
		integrator_var_to_inertial_pos();
	}
	integrator_update_acceleration();
	integrator_to_jacobi_acc();
	if (N_megno){
		integrator_var_to_jacobi_acc();
	}
	integrator_interaction(b);
	integrator_kepler_drift(a);
}

static void integrator_apply_corrector(double inv){
	if (integrator_whfast_corrector==3){
		// Third order corrector
		integrator_whfast_corrector_Z(a_1*dt,-inv*b_31*dt);
		integrator_whfast_corrector_Z(-a_1*dt,inv*b_31*dt);
	}
	if (integrator_whfast_corrector==5){
		// Fifth order corrector
		integrator_whfast_corrector_Z(-a_2*dt,-inv*b_51*dt);
		integrator_whfast_corrector_Z(-a_1*dt,-inv*b_52*dt);
		integrator_whfast_corrector_Z(a_1*dt,inv*b_52*dt);
		integrator_whfast_corrector_Z(a_2*dt,inv*b_51*dt);
	}
	if (integrator_whfast_corrector==7){
		// Seventh order corrector
		integrator_whfast_corrector_Z(-a_3*dt,-inv*b_71*dt);
		integrator_whfast_corrector_Z(-a_2*dt,-inv*b_72*dt);
		integrator_whfast_corrector_Z(-a_1*dt,-inv*b_73*dt);
		integrator_whfast_corrector_Z(a_1*dt,inv*b_73*dt);
		integrator_whfast_corrector_Z(a_2*dt,inv*b_72*dt);
		integrator_whfast_corrector_Z(a_3*dt,inv*b_71*dt);
	}
	if (integrator_whfast_corrector==11){
		// Eleventh order corrector
		integrator_whfast_corrector_Z(-a_5*dt,-inv*b_111*dt);
		integrator_whfast_corrector_Z(-a_4*dt,-inv*b_112*dt);
		integrator_whfast_corrector_Z(-a_3*dt,-inv*b_113*dt);
		integrator_whfast_corrector_Z(-a_2*dt,-inv*b_114*dt);
		integrator_whfast_corrector_Z(-a_1*dt,-inv*b_115*dt);
		integrator_whfast_corrector_Z(a_1*dt,inv*b_115*dt);
		integrator_whfast_corrector_Z(a_2*dt,inv*b_114*dt);
		integrator_whfast_corrector_Z(a_3*dt,inv*b_113*dt);
		integrator_whfast_corrector_Z(a_4*dt,inv*b_112*dt);
		integrator_whfast_corrector_Z(a_5*dt,inv*b_111*dt);
	}
}

void integrator_whfast_part1(){
	gravity_ignore_10 = 1;
	if (integrator_allocated_N != N){
		integrator_allocated_N = N;
		p_j = realloc(p_j,sizeof(struct particle)*N);
		eta = realloc(eta,sizeof(double)*(N-N_megno));
		etai= realloc(etai,sizeof(double)*(N-N_megno));
		integrator_whfast_recalculate_jacobi_this_timestep = 1;
	}
	// Only recalculate Jacobi coordinates if needed
	if (integrator_whfast_safe_mode || integrator_whfast_recalculate_jacobi_this_timestep){
		if (integrator_whfast_is_synchronized==0){
			integrator_whfast_synchronize();
			if (integrator_whfast_recalculate_jacobi_but_not_synchronized_warning==0){
				fprintf(stderr,"\n\033[1mWarning!\033[0m Recalculating Jacobi coordinates but pos/vel were not synchronized before.\n");
				integrator_whfast_recalculate_jacobi_but_not_synchronized_warning++;
			}
		}
		eta[0] = particles[0].m;
		etai[0] = 1./eta[0];
		p_j[0].m = particles[0].m;
		for (unsigned int i=1;i<N-N_megno;i++){
			eta[i] = eta[i-1] + particles[i].m;
			etai[i] = 1./eta[i];
			p_j[i].m = particles[i].m;
		}
		for (unsigned int i=N_megno;i<N;i++){
			p_j[i].m = particles[i].m;
		}
		Mtotal  = eta[N-N_megno-1];
		Mtotali = etai[N-N_megno-1];
		integrator_whfast_recalculate_jacobi_this_timestep = 0;
		integrator_to_jacobi_posvel();
		if (N_megno){
			integrator_var_to_jacobi_posvel();
		}
	}
	double _dt2 = dt/2.;
	if (integrator_whfast_is_synchronized){
		// First half DRIFT step
		if (integrator_whfast_corrector){
			integrator_apply_corrector(1.);
		}
		integrator_kepler_drift(_dt2);	// half timestep
	}else{
		// Combined DRIFT step
		integrator_kepler_drift(dt);	// full timestep
	}
	// Prepare coordinates for KICK step
	if (integrator_force_is_velocitydependent){
		integrator_to_inertial_posvel();
	}else{
		integrator_to_inertial_pos();
	}
	
	if (N_megno){
		p_j[N_megno].x += _dt2*p_j[N_megno].vx;
		p_j[N_megno].y += _dt2*p_j[N_megno].vy;
		p_j[N_megno].z += _dt2*p_j[N_megno].vz;
		if (integrator_force_is_velocitydependent){
			integrator_var_to_inertial_posvel();
		}else{
			integrator_var_to_inertial_pos();
		}
	}

	t+=_dt2;
}

void integrator_whfast_synchronize(void){
	if (integrator_whfast_is_synchronized == 0){
		integrator_kepler_drift(dt/2.);
		if (integrator_whfast_corrector){
			integrator_apply_corrector(-1.);
		}
		integrator_to_inertial_posvel();
		integrator_whfast_is_synchronized = 1;
	}
}

void integrator_whfast_part2(void){
	integrator_to_jacobi_acc();
	if (N_megno){
		integrator_var_to_jacobi_acc();
	}
	integrator_interaction(dt);

	double _dt2 = dt/2.;
	integrator_whfast_is_synchronized = 0;
	if (integrator_whfast_safe_mode || N_megno){
		integrator_whfast_synchronize();
	}
	
	t+=_dt2;

	if (N_megno){
		p_j[N_megno].x += _dt2*p_j[N_megno].vx;
		p_j[N_megno].y += _dt2*p_j[N_megno].vy;
		p_j[N_megno].z += _dt2*p_j[N_megno].vz;
		integrator_var_to_inertial_posvel();
		gravity_calculate_variational_acceleration();
		// Add additional acceleration term for MEGNO calculation
		int i = N-N_megno;
		int j = N-N_megno+1;
		const double dx = particles[i-N/2].x - particles[j-N/2].x;
		const double dy = particles[i-N/2].y - particles[j-N/2].y;
		const double dz = particles[i-N/2].z - particles[j-N/2].z;
		const double r2 = dx*dx + dy*dy + dz*dz + softening*softening;
		const double r  = sqrt(r2);
		const double r3inv = 1./(r2*r);
		const double r5inv = 3.*r3inv/r2;
		const double ddx = particles[i].x - particles[j].x;
		const double ddy = particles[i].y - particles[j].y;
		const double ddz = particles[i].z - particles[j].z;
		const double Gmi = G * particles[i].m;
		const double Gmj = G * particles[j].m;
		const double dax =   ddx * ( dx*dx*r5inv - r3inv )
				   + ddy * ( dx*dy*r5inv )
				   + ddz * ( dx*dz*r5inv );
		const double day =   ddx * ( dy*dx*r5inv )
				   + ddy * ( dy*dy*r5inv - r3inv )
				   + ddz * ( dy*dz*r5inv );
		const double daz =   ddx * ( dz*dx*r5inv )
				   + ddy * ( dz*dy*r5inv )
				   + ddz * ( dz*dz*r5inv - r3inv );
		
		particles[i].ax += Gmj * dax;
		particles[i].ay += Gmj * day;
		particles[i].az += Gmj * daz;
		
		particles[j].ax -= Gmi * dax;
		particles[j].ay -= Gmi * day;
		particles[j].az -= Gmi * daz;

		// Update MEGNO in middle of timestep as we need synchonized x/v/a.
		double dY = dt * 2. * t * tools_megno_deltad_delta();
		tools_megno_update(dY);
	}
}
	
void integrator_whfast_reset(void){
	integrator_whfast_corrector = 0;
	integrator_whfast_is_synchronized = 1;
	integrator_whfast_safe_mode = 1;
	integrator_whfast_recalculate_jacobi_this_timestep = 0;
	integrator_allocated_N = 0;
	integrator_timestep_warning = 0;
	integrator_whfast_recalculate_jacobi_but_not_synchronized_warning = 0;
	free(p_j);
	p_j = NULL;
	free(eta);
	eta = NULL;
	free(etai);
	etai = NULL;
}
