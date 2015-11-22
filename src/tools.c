/**
 * @file 	tools.c
 * @brief 	Tools for creating distributions.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
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
#include <sys/time.h>
#ifndef LIBREBOUNDX
#include "particle.h"
#endif // LIBREBOUNDX
#include "rebound.h"
#include "tools.h"


void reb_tools_init_srand(void){
	struct timeval tim;
	gettimeofday(&tim, NULL);
	srand ( tim.tv_usec + getpid());
}

double reb_random_uniform(double min, double max){
	return ((double)rand())/((double)(RAND_MAX))*(max-min)+min;
}


double reb_random_powerlaw(double min, double max, double slope){
	double y = reb_random_uniform(0., 1.);
	return pow( (pow(max,slope+1.)-pow(min,slope+1.))*y+pow(min,slope+1.), 1./(slope+1.));
}

double reb_random_normal(double variance){
	double v1,v2,rsq=1.;
	while(rsq>=1. || rsq<1.0e-12){
		v1=2.*((double)rand())/((double)(RAND_MAX))-1.0;
		v2=2.*((double)rand())/((double)(RAND_MAX))-1.0;
		rsq=v1*v1+v2*v2;
	}
	// Note: This gives another random variable for free, but we'll throw it away for simplicity and for thread-safety.
	return 	v1*sqrt(-2.*log(rsq)/rsq*variance);
}

double reb_random_rayleigh(double sigma){
	double y = reb_random_uniform(0.,1.);
	return sigma*sqrt(-2*log(y));
}

/// Other helper routines
double reb_tools_energy(struct reb_simulation* r){
	const int N = r->N;
	const struct reb_particle* restrict const particles = r->particles;
	const int N_var = r->N_var;
	double e_kin = 0.;
	double e_pot = 0.;
	for (int i=0;i<N-N_var;i++){
		struct reb_particle pi = particles[i];
		e_kin += 0.5 * pi.m * (pi.vx*pi.vx + pi.vy*pi.vy + pi.vz*pi.vz);
		for (int j=i+1;j<N-N_var;j++){
			struct reb_particle pj = particles[j];
			double dx = pi.x - pj.x;
			double dy = pi.y - pj.y;
			double dz = pi.z - pj.z;
			e_pot -= r->G*pj.m*pi.m/sqrt(dx*dx + dy*dy + dz*dz);
		}
	}
	return e_kin +e_pot;
}

void reb_move_to_com(struct reb_simulation* const r){
	const int N = r->N;
	struct reb_particle* restrict const particles = r->particles;
	struct reb_particle com = reb_get_com(r);
	for (int i=0;i<N;i++){
		particles[i].x  -= com.x;
		particles[i].y  -= com.y;
		particles[i].z  -= com.z;
		particles[i].vx -= com.vx;
		particles[i].vy -= com.vy;
		particles[i].vz -= com.vz;
	}
}

struct reb_particle reb_get_com(struct reb_simulation* r){
	struct reb_particle com = {.m=0, .x=0, .y=0, .z=0, .vx=0, .vy=0, .vz=0};
	const int N = r->N;
	struct reb_particle* restrict const particles = r->particles;
	for (int i=0;i<N;i++){
		com = reb_get_com_of_pair(com, particles[i]);
	}
	return com;
}

struct reb_particle reb_get_com_of_pair(struct reb_particle p1, struct reb_particle p2){
	p1.x   = p1.x*p1.m + p2.x*p2.m;		
	p1.y   = p1.y*p1.m + p2.y*p2.m;
	p1.z   = p1.z*p1.m + p2.z*p2.m;
	p1.vx  = p1.vx*p1.m + p2.vx*p2.m;
	p1.vy  = p1.vy*p1.m + p2.vy*p2.m;
	p1.vz  = p1.vz*p1.m + p2.vz*p2.m;
	p1.m  += p2.m;
	if (p1.m>0.){
		p1.x  /= p1.m;
		p1.y  /= p1.m;
		p1.z  /= p1.m;
		p1.vx /= p1.m;
		p1.vy /= p1.m;
		p1.vz /= p1.m;
	}
	return p1;
}

int reb_get_particle_index(struct reb_particle* p){
	struct reb_simulation* r = p->sim;
	int i = 0;
	int N = r->N-r->N_var;
	while(&r->particles[i] != p){
		i++;
		if(i>=N){
			return -1;	// p not in simulation.  Shouldn't happen unless you mess with p.sim after creating the particle
		}	
	}
	return i;
}

struct reb_particle reb_get_jacobi_com(struct reb_particle* p){
	int p_index = reb_get_particle_index(p);
	struct reb_simulation* r = p->sim;
	struct reb_particle com = r->particles[0];
	for(int i=1; i<p_index; i++){
		com = reb_get_com_of_pair(com, r->particles[i]);
	}
	return com;
}
	
#ifndef LIBREBOUNDX
void reb_tools_init_plummer(struct reb_simulation* r, int _N, double M, double R) {
	// Algorithm from:	
	// http://adsabs.harvard.edu/abs/1974A%26A....37..183A
	
	double E = 3./64.*M_PI*M*M/R;
	for (int i=0;i<_N;i++){
		struct reb_particle star = {0};
		double _r = pow(pow(reb_random_uniform(0,1),-2./3.)-1.,-1./2.);
		double x2 = reb_random_uniform(0,1);
		double x3 = reb_random_uniform(0,2.*M_PI);
		star.z = (1.-2.*x2)*_r;
		star.x = sqrt(_r*_r-star.z*star.z)*cos(x3);
		star.y = sqrt(_r*_r-star.z*star.z)*sin(x3);
		double x5,g,q;
		do{
			x5 = reb_random_uniform(0.,1.);
			q = reb_random_uniform(0.,1.);
			g = q*q*pow(1.-q*q,7./2.);
		}while(0.1*x5>g);
		double ve = pow(2.,1./2.)*pow(1.+_r*_r,-1./4.);
		double v = q*ve;
		double x6 = reb_random_uniform(0.,1.);
		double x7 = reb_random_uniform(0.,2.*M_PI);
		star.vz = (1.-2.*x6)*v;
		star.vx = sqrt(v*v-star.vz*star.vz)*cos(x7);
		star.vy = sqrt(v*v-star.vz*star.vz)*sin(x7);
		
		star.x *= 3.*M_PI/64.*M*M/E;
		star.y *= 3.*M_PI/64.*M*M/E;
		star.z *= 3.*M_PI/64.*M*M/E;
		
		star.vx *= sqrt(E*64./3./M_PI/M);
		star.vy *= sqrt(E*64./3./M_PI/M);
		star.vz *= sqrt(E*64./3./M_PI/M);

		star.m = M/(double)_N;

		reb_add(r, star);
	}
}
#endif // LIBREBOUNDX

double mod2pi(double f){
	while(f < 0.){
		f += 2*M_PI;
	}
	while(f > 0.){
		f -= 2*M_PI;
	}
	return f;
}

double reb_M_to_E(double e, double M){
	double E;
	if(e < 1.){
		E = e < 0.8 ? M : M_PI;
		double F = E - e*sin(M) - M;
		for(int i=0; i<100; i++){
			E = E - F/(1.-e*cos(E));
			F = E - e*sin(E) - M;
			if(fabs(F) < 1.e-16){
				break;
			}
		}
		E = mod2pi(E);
		return E;
	}
	else{
		E = M;

		double F = E - e*sinh(E) - M;
		for(int i=0; i<100; i++){
			E = E - F/(1.0 - e*cosh(E));
			F = E - e*sinh(E) - M;
			if(fabs(F) < 1.e-16){
				break;
			}
		}
		E = mod2pi(E);
		return E;
	}
}

double reb_tools_M_to_f(double e, double M){

	double E = reb_M_to_E(e, M);
	if(e > 1.){
		return 2.*atan(sqrt((1.+e)/(e-1.))*tanh(0.5*E));
	}
	else{
		return 2*atan(sqrt((1.+e)/(1.-e))*tan(0.5*E));
	}
}

struct reb_particle reb_tools_orbit2d_to_particle(double G, struct reb_particle primary, double m, double a, double e, double omega, double f){
	double Omega = 0.;
	double inc = 0.;
	return reb_tools_orbit_to_particle(G, primary, m, a, e, inc, Omega, omega, f);
}

static const struct reb_particle reb_particle_nan = {.x = NAN, .y = NAN, .z = NAN, .vx = NAN, .vy = NAN, .vz = NAN, .ax = NAN, .ay = NAN, .az = NAN, .m = NAN, .r = NAN, .lastcollision = NAN, .c = NULL, .id = -1, .ap = NULL, .sim = NULL};

struct reb_particle reb_tools_orbit_to_particle_err(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f, int* err){
	if(e == 1.){
		*err = 1; 		// Can't initialize a radial orbit with orbital elements.
		return reb_particle_nan;
	}
	if(e < 0.){
		*err = 2; 		// Eccentricity must be greater than or equal to zero.
		return reb_particle_nan;
	}
	if(e > 1.){
		if(a > 0.){
			*err = 3; 	// Bound orbit (a > 0) must have e < 1. 
			return reb_particle_nan;
		}
	}
	else{
		if(a < 0.){
			*err =4; 	// Unbound orbit (a < 0) must have e > 1.
			return reb_particle_nan;
		}
	}
	if(e*cos(f) < -1.){
		*err = 5;		// Unbound orbit can't have f set beyond the range allowed by the asymptotes set by the parabola.
		return reb_particle_nan;
	}

	struct reb_particle p = {0};
	p.m = m;
	double r = a*(1-e*e)/(1 + e*cos(f));
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); // in this form it works for elliptical and hyperbolic orbits

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	// Murray & Dermott Eq 2.122
	p.x = primary.x + r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
	p.y = primary.y + r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
	p.z = primary.z + r*(so*cf+co*sf)*si;

	// Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
	p.vx = primary.vx + v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
	p.vy = primary.vy + v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
	p.vz = primary.vz + v0*((e+cf)*co*si - sf*si*so);
	
	p.ax = 0; 	p.ay = 0; 	p.az = 0;

	return p;
}

struct reb_particle reb_tools_orbit_to_particle(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	int err;
	return reb_tools_orbit_to_particle_err(G, primary, m, a, e, inc, Omega, omega, f, &err);
}

static const struct reb_orbit reb_orbit_nan = {.r = NAN, .v = NAN, .h = NAN, .P = NAN, .n = NAN, .a = NAN, .e = NAN, .inc = NAN, .Omega = NAN, .omega = NAN, .pomega = NAN, .f = NAN, .M = NAN, .l = NAN};

#define MIN_REL_ERROR 1.0e-12	///< Close to smallest relative floating point number, used for orbit calculation
#define TINY 1.E-308 		///< Close to smallest representable floating point number, used for orbit calculation
#define MIN_INC 1.e-8		///< Below this inclination, the broken angles pomega and theta equal the corresponding 
							///< unbroken angles to within machine precision, so a practical boundary for planar orbits
							//
// returns acos(num/denom), using disambiguator to tell which quadrant to return.  
// will return 0 or pi appropriately if num is larger than denom by machine precision
// and will return 0 if denom is exactly 0.

static double acos2(double num, double denom, double disambiguator){
	double val;
	double cosine = num/denom;
	if(cosine > -1. && cosine < 1.){
		val = acos(cosine);
		if(disambiguator < 0.){
			val = - val;
		}
	}
	else{
		val = (cosine <= -1.) ? M_PI : 0.;
	}
	return val;
}

struct reb_orbit reb_tools_particle_to_orbit_err(double G, struct reb_particle p, struct reb_particle primary, int* err){
	struct reb_orbit o;
	if (primary.m <= TINY){	
		*err = 1;			// primary has no mass.
		return reb_orbit_nan;
	}
	double mu,dx,dy,dz,dvx,dvy,dvz,vsquared,vcircsquared,vdiffsquared;
	double hx,hy,hz,vr,rvr,muinv,ex,ey,ez,nx,ny,n,ea;
	mu = G*(p.m+primary.m);
	dx = p.x - primary.x;
	dy = p.y - primary.y;
	dz = p.z - primary.z;
	dvx = p.vx - primary.vx;
	dvy = p.vy - primary.vy;
	dvz = p.vz - primary.vz;
	o.r = sqrt ( dx*dx + dy*dy + dz*dz );
	
	vsquared = dvx*dvx + dvy*dvy + dvz*dvz;
	o.v = sqrt(vsquared);
	vcircsquared = mu/o.r;	
	o.a = -mu/( vsquared - 2.*vcircsquared );	// semi major axis
	
	hx = (dy*dvz - dz*dvy); 					//angular momentum vector
	hy = (dz*dvx - dx*dvz);
	hz = (dx*dvy - dy*dvx);
	o.h = sqrt ( hx*hx + hy*hy + hz*hz );		// abs value of angular momentum

	vdiffsquared = vsquared - vcircsquared;	
	if(o.r <= TINY){							
		*err = 2;									// particle is on top of primary
		return reb_orbit_nan;
	}
	vr = (dx*dvx + dy*dvy + dz*dvz)/o.r;	
	rvr = o.r*vr;
	muinv = 1./mu;

	ex = muinv*( vdiffsquared*dx - rvr*dvx );
	ey = muinv*( vdiffsquared*dy - rvr*dvy );
	ez = muinv*( vdiffsquared*dz - rvr*dvz );
 	o.e = sqrt( ex*ex + ey*ey + ez*ez );		// eccentricity
	o.n = o.a/fabs(o.a)*sqrt(fabs(mu/(o.a*o.a*o.a)));	// mean motion (negative if hyperbolic)
	o.P = 2*M_PI/o.n;									// period (negative if hyperbolic)

	o.inc = acos2(hz, o.h, 1.);			// cosi = dot product of h and z unit vectors.  Always in [0,pi], so pass dummy disambiguator
										// will = 0 if h is 0.

	nx = -hy;							// vector pointing along the ascending node = zhat cross h
	ny =  hx;		
	n = sqrt( nx*nx + ny*ny );

	// Omega, pomega and theta are measured from x axis, so we can always use y component to disambiguate if in the range [0,pi] or [pi,2pi]
	o.Omega = acos2(nx, n, ny);			// cos Omega is dot product of x and n unit vectors. Will = 0 if i=0.
	
	ea = acos2(1.-o.r/o.a, o.e, vr);	// from definition of eccentric anomaly.  If vr < 0, must be going from apo to peri, so ea = [pi, 2pi] so ea = -acos(cosea)
	o.M = ea - o.e*sin(ea);						// mean anomaly (Kepler's equation)

	// in the near-planar case, the true longitude is always well defined for the position, and pomega for the pericenter if e!= 0
	// we therefore calculate those and calculate the remaining angles from them
	if(o.inc < MIN_INC || o.inc > M_PI - MIN_INC){	// nearly planar.  Use longitudes rather than angles referenced to node for numerical stability.
		o.pomega = acos2(ex, o.e, ey);		// cos pomega is dot product of x and e unit vectors.  Will = 0 if e=0.
		o.theta = acos2(dx, o.r, dy);			// cos theta is dot product of x and r vectors (true longitude).  Will = 0 if e = 0.
		if(o.inc < M_PI/2.){
			o.omega = o.pomega - o.Omega;
			o.f = o.theta - o.pomega;
			o.l = o.pomega + o.M;
		}
		else{
			o.omega = o.Omega - o.pomega;
			o.f = o.pomega - o.theta;
			o.l = o.pomega - o.M;
		}
	}
	// in the non-planar case, we can't calculate the broken angles from vectors like above.  omega+f is always well defined, and omega if e!=0
	else{
		double wpf = acos2(nx*dx + ny*dy, n*o.r, dz);	// omega plus f.  Both angles measured in orbital plane, and always well defined for i!=0.
		o.omega = acos2(nx*ex + ny*ey, n*o.e, ez);
		if(o.inc < M_PI/2.){
			o.pomega = o.Omega + o.omega;
			o.f = wpf - o.omega;
			o.theta = o.Omega + wpf;
			o.l = o.pomega + o.M;
		}
		else{
			o.pomega = o.Omega - o.omega;
			o.f = wpf - o.omega;
			o.theta = o.Omega - wpf;
			o.l = o.pomega - o.M;
		}
	}

	return o;
}


struct reb_orbit reb_tools_particle_to_orbit(double G, struct reb_particle p, struct reb_particle primary){
	int err;
	return reb_tools_particle_to_orbit_err(G, p, primary, &err);
}

/**************************
 * MEGNO Routines         */

#ifndef LIBREBOUNDX
void reb_tools_megno_init(struct reb_simulation* const r, double delta){
	int N_var = r->N;
	r->calculate_megno = 1;
	r->megno_Ys = 0.;
	r->megno_Yss = 0.;
	r->megno_cov_Yt = 0.;
	r->megno_var_t = 0.;
	r->megno_n = 0;
	r->megno_mean_Y = 0;
	r->megno_mean_t = 0;
        for (int i=0;i<N_var;i++){ 
                struct reb_particle megno = {
			.m  = r->particles[i].m,
			.x  = reb_random_normal(1.),
			.y  = reb_random_normal(1.),
			.z  = reb_random_normal(1.),
			.vx = reb_random_normal(1.),
			.vy = reb_random_normal(1.),
			.vz = reb_random_normal(1.) };
		double deltad = delta/sqrt(megno.x*megno.x + megno.y*megno.y + megno.z*megno.z + megno.vx*megno.vx + megno.vy*megno.vy + megno.vz*megno.vz); // rescale
		megno.x *= deltad;
		megno.y *= deltad;
		megno.z *= deltad;
		megno.vx *= deltad;
		megno.vy *= deltad;
		megno.vz *= deltad;

                reb_add(r, megno);
        }
	r->N_var = N_var;
}
double reb_tools_calculate_megno(struct reb_simulation* r){ // Returns the MEGNO <Y>
	if (r->t==0.) return 0.;
	return r->megno_Yss/r->t;
}
double reb_tools_calculate_lyapunov(struct reb_simulation* r){ // Returns the largest Lyapunov characteristic number (LCN), or maximal Lyapunov exponent
	if (r->t==0.) return 0.;
	return r->megno_cov_Yt/r->megno_var_t;
}
double reb_tools_megno_deltad_delta(struct reb_simulation* const r){
	const struct reb_particle* restrict const particles = r->particles;
	const int N = r->N;
	const int N_var = r->N_var;
        double deltad = 0;
        double delta2 = 0;
        for (int i=N-N_var;i<N;i++){
                deltad += particles[i].vx * particles[i].x; 
                deltad += particles[i].vy * particles[i].y; 
                deltad += particles[i].vz * particles[i].z; 
                deltad += particles[i].ax * particles[i].vx; 
                deltad += particles[i].ay * particles[i].vy; 
                deltad += particles[i].az * particles[i].vz; 
                delta2 += particles[i].x  * particles[i].x; 
                delta2 += particles[i].y  * particles[i].y;
                delta2 += particles[i].z  * particles[i].z;
                delta2 += particles[i].vx * particles[i].vx; 
                delta2 += particles[i].vy * particles[i].vy;
                delta2 += particles[i].vz * particles[i].vz;
        }
        return deltad/delta2;
}

void reb_tools_megno_update(struct reb_simulation* r, double dY){
	// Calculate running Y(t)
	r->megno_Ys += dY;
	double Y = r->megno_Ys/r->t;
	// Calculate averge <Y> 
	r->megno_Yss += Y * r->dt;
	// Update covariance of (Y,t) and variance of t
	r->megno_n++;
	double _d_t = r->t - r->megno_mean_t;
	r->megno_mean_t += _d_t/(double)r->megno_n;
	double _d_Y = reb_tools_calculate_megno(r) - r->megno_mean_Y;
	r->megno_mean_Y += _d_Y/(double)r->megno_n;
	r->megno_cov_Yt += ((double)r->megno_n-1.)/(double)r->megno_n 
					*(r->t-r->megno_mean_t)
					*(reb_tools_calculate_megno(r)-r->megno_mean_Y);
	r->megno_var_t  += ((double)r->megno_n-1.)/(double)r->megno_n 
					*(r->t-r->megno_mean_t)
					*(r->t-r->megno_mean_t);
}
#endif // LIBREBOUNDX
