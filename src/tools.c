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
#include "particle.h"
#include "rebound.h"
#include "tools.h"


void reb_tools_init_srand(){
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

struct reb_particle reb_tools_init_orbit2d(double G, double M, double m, double a, double e, double omega, double f){
	struct reb_particle p = {0};
	p.m = m;
	double r = a*(1.-e*e)/(1.+e*cos(f));
	double n = sqrt(G*(m+M)/(a*a*a));
	double tx = r*cos(f);
	double ty = r*sin(f);
	p.z = 0;
	double tvx = -n*a/sqrt(1.-e*e)*sin(f);
	double tvy =  n*a/sqrt(1.-e*e)*(e+cos(f));
	p.vz = 0;
	p.ax = 0; p.ay = 0; p.az = 0;
	
	p.x  =  cos(omega)*tx  - sin(omega)*ty;
	p.y  =  sin(omega)*tx  + cos(omega)*ty;
	
	p.vx =  cos(omega)*tvx - sin(omega)*tvy;
	p.vy =  sin(omega)*tvx + cos(omega)*tvy;

	return p;
}

struct reb_particle reb_tools_init_orbit3d(double G, double M, double m, double a, double e, double i, double Omega, double omega, double f){
	struct reb_particle p = {0};
	p.m = m;
	double r = a*(1-e*e)/(1 + e*cos(f));

	// Murray & Dermott Eq 2.122
	p.x  = r*(cos(Omega)*cos(omega+f) - sin(Omega)*sin(omega+f)*cos(i));
	p.y  = r*(sin(Omega)*cos(omega+f) + cos(Omega)*sin(omega+f)*cos(i));
	p.z  = r*sin(omega+f)*sin(i);

	double n = sqrt(G*(m+M)/(a*a*a));

	// Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
	p.vx = (n*a/sqrt(1-e*e))*((e+cos(f))*(-cos(i)*cos(omega)*sin(Omega) - cos(Omega)*sin(omega)) - sin(f)*(cos(omega)*cos(Omega) - cos(i)*sin(omega)*sin(Omega)));
	p.vy = (n*a/sqrt(1-e*e))*((e+cos(f))*(cos(i)*cos(omega)*cos(Omega) - sin(omega)*sin(Omega)) - sin(f)*(cos(omega)*sin(Omega) + cos(i)*cos(Omega)*sin(omega)));
	p.vz = (n*a/sqrt(1-e*e))*((e+cos(f))*cos(omega)*sin(i) - sin(f)*sin(i)*sin(omega));

	p.ax = 0; 	p.ay = 0; 	p.az = 0;

	return p;
}

static const struct reb_orbit reb_orbit_nan = {.a = NAN, .r = NAN, .h = NAN, .P = NAN, .l = NAN, .e = NAN, .inc = NAN, .Omega = NAN, .omega = NAN, .f = NAN};

#define MIN_REL_ERROR 1.0e-12	///< Close to smallest relative floating point number, used for orbit calculation
#define TINY 1.E-308 		///< Close to smallest representable floating point number, used for orbit calculation

struct reb_orbit reb_tools_p2orbit(double G, struct reb_particle p, struct reb_particle primary){
	struct reb_orbit o;
	if (primary.m <= TINY){	
		return reb_orbit_nan;
	}
	double h0,h1,h2,e0,e1,e2,n0,n1,n,er,vr,mu,ea,dx,dy,dz,dvx,dvy,dvz,v,cosf,cosea;
	mu = G*(p.m+primary.m);
	dx = p.x - primary.x;
	dy = p.y - primary.y;
	dz = p.z - primary.z;
	dvx = p.vx - primary.vx;
	dvy = p.vy - primary.vy;
	dvz = p.vz - primary.vz;
	h0 = (dy*dvz - dz*dvy); 			//angular momentum vector
	h1 = (dz*dvx - dx*dvz);
	h2 = (dx*dvy - dy*dvx);
	o.h = sqrt ( h0*h0 + h1*h1 + h2*h2 );		// abs value of angular moment 
	v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
	o.r = sqrt ( dx*dx + dy*dy + dz*dz );
	if(o.r <= TINY){
		return reb_orbit_nan;
	}
	if (o.h/(o.r*v) <= MIN_REL_ERROR){
		return reb_orbit_nan;
	}
	vr = (dx*dvx + dy*dvy + dz*dvz)/o.r;
	e0 = 1./mu*( (v*v-mu/o.r)*dx - o.r*vr*dvx );
	e1 = 1./mu*( (v*v-mu/o.r)*dy - o.r*vr*dvy );
	e2 = 1./mu*( (v*v-mu/o.r)*dz - o.r*vr*dvz );
 	o.e = sqrt( e0*e0 + e1*e1 + e2*e2 );		// eccentricity
	o.a = -mu/( v*v - 2.*mu/o.r );			// semi major axis

	o.P = o.a/fabs(o.a)*2.*M_PI*sqrt(fabs(o.a*o.a*o.a/mu));		// period (negative if hyperbolic)
	o.inc = acos( h2/o.h ) ;				// inclination (wrt xy-plane)   -  Note if pi/2 < i < pi then the orbit is retrograde
	n0 = -h1;					// vector of nodes lies in xy plane => no z component
	n1 =  h0;		
	n = sqrt( n0*n0 + n1*n1 );
	er = dx*e0 + dy*e1 + dz*e2;
	if (n/(o.r*v) <= MIN_REL_ERROR || o.inc <= MIN_REL_ERROR){			// we are in the xy plane
		o.Omega=0.;
		if (o.e <= MIN_REL_ERROR){              // omega not defined for circular orbit
			o.omega = 0.;
		}
		else{
			if (e1>=0.){
				o.omega=acos(e0/o.e);
			}
			else{
				o.omega = 2.*M_PI-acos(e0/o.e);
			}
		}
	}
	else{
		if (o.e <= MIN_REL_ERROR){
			o.omega = 0.;
		}
		else{
			if (e2>=0.){                        // omega=0 if perictr at asc node
				o.omega=acos(( n0*e0 + n1*e1 )/(n*o.e));
			}
			else{
				o.omega=2.*M_PI-acos(( n0*e0 + n1*e1 )/(n*o.e));
			}
		}

		if (n1>=0.){
			o.Omega = acos(n0/n);
		}
		else{
			o.Omega=2.*M_PI-acos(n0/n); // Omega=longitude of asc node
		}								// taken in xy plane from x axis
	}	
	if (o.e<=MIN_REL_ERROR){            // circular orbit
		o.f=0.;                         // f has no meaning
		o.l=0.;
	}
	else{
		cosf = er/(o.e*o.r);
		cosea = (1.-o.r/o.a)/o.e;
		
		if (-1.<=cosf && cosf<=1.){     // failsafe
			o.f = acos(cosf);
		}
		else{
			o.f = M_PI/2.*(1.-cosf);
		}
		
		if (-1.<=cosea && cosea<=1.){
			ea  = acos(cosea);
		}
		else{
			ea = M_PI/2.*(1.-cosea);
		}
		
		if (vr<0.){
			o.f=2.*M_PI-o.f;
			ea =2.*M_PI-ea;
		}
		
		o.l = ea -o.e*sin(ea) + o.omega+ o.Omega;  // mean longitude
	}
	return o;
}

/**************************
 * MEGNO Routines         */

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
