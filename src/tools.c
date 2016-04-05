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
	if(slope == -1) return exp(y*log(max/min) + log(min));
    else return pow( (pow(max,slope+1.)-pow(min,slope+1.))*y+pow(min,slope+1.), 1./(slope+1.));
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
double reb_tools_energy(const struct reb_simulation* const r){
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
    const int N_real = r->N - r->N_var;
	struct reb_particle* restrict const particles = r->particles;
	struct reb_particle com = reb_get_com(r);
    // First do second order
    for (int v=0;v<r->var_config_N;v++){
        int index = r->var_config[v].index;
        if (r->var_config[v].testparticle>=0){
            // Test particles do not affect the COM
        }else{
            if (r->var_config[v].order==2){
                struct reb_particle com_shift = {0};
                int index_1st_order_a = r->var_config[v].index_1st_order_a;
                int index_1st_order_b = r->var_config[v].index_1st_order_b;
                double dma = 0.;
                double dmb = 0.;
                double ddm = 0.;
                for (int i=0;i<N_real;i++){
                    dma += particles[i+index_1st_order_a].m;
                    dmb += particles[i+index_1st_order_b].m;
                    ddm += particles[i+index].m;
                }
                for (int i=0;i<N_real;i++){
                    com_shift.x  += particles[i+index].x /com.m * particles[i].m; 
                    com_shift.y  += particles[i+index].y /com.m * particles[i].m; 
                    com_shift.z  += particles[i+index].z /com.m * particles[i].m; 
                    com_shift.vx += particles[i+index].vx/com.m * particles[i].m; 
                    com_shift.vy += particles[i+index].vy/com.m * particles[i].m; 
                    com_shift.vz += particles[i+index].vz/com.m * particles[i].m; 
                    
                    com_shift.x  += particles[i+index_1st_order_a].x  /com.m * particles[i+index_1st_order_b].m; 
                    com_shift.y  += particles[i+index_1st_order_a].y  /com.m * particles[i+index_1st_order_b].m; 
                    com_shift.z  += particles[i+index_1st_order_a].z  /com.m * particles[i+index_1st_order_b].m; 
                    com_shift.vx += particles[i+index_1st_order_a].vx /com.m * particles[i+index_1st_order_b].m; 
                    com_shift.vy += particles[i+index_1st_order_a].vy /com.m * particles[i+index_1st_order_b].m; 
                    com_shift.vz += particles[i+index_1st_order_a].vz /com.m * particles[i+index_1st_order_b].m; 
                    
                    com_shift.x  -= particles[i+index_1st_order_a].x  * particles[i].m/com.m/com.m*dmb; 
                    com_shift.y  -= particles[i+index_1st_order_a].y  * particles[i].m/com.m/com.m*dmb; 
                    com_shift.z  -= particles[i+index_1st_order_a].z  * particles[i].m/com.m/com.m*dmb; 
                    com_shift.vx -= particles[i+index_1st_order_a].vx * particles[i].m/com.m/com.m*dmb; 
                    com_shift.vy -= particles[i+index_1st_order_a].vy * particles[i].m/com.m/com.m*dmb; 
                    com_shift.vz -= particles[i+index_1st_order_a].vz * particles[i].m/com.m/com.m*dmb; 
                    
                    com_shift.x  += particles[i+index_1st_order_b].x  /com.m * particles[i+index_1st_order_a].m; 
                    com_shift.y  += particles[i+index_1st_order_b].y  /com.m * particles[i+index_1st_order_a].m; 
                    com_shift.z  += particles[i+index_1st_order_b].z  /com.m * particles[i+index_1st_order_a].m; 
                    com_shift.vx += particles[i+index_1st_order_b].vx /com.m * particles[i+index_1st_order_a].m; 
                    com_shift.vy += particles[i+index_1st_order_b].vy /com.m * particles[i+index_1st_order_a].m; 
                    com_shift.vz += particles[i+index_1st_order_b].vz /com.m * particles[i+index_1st_order_a].m; 
                   
                    com_shift.x  += particles[i].x  /com.m * particles[i+index].m; 
                    com_shift.y  += particles[i].y  /com.m * particles[i+index].m; 
                    com_shift.z  += particles[i].z  /com.m * particles[i+index].m; 
                    com_shift.vx += particles[i].vx /com.m * particles[i+index].m; 
                    com_shift.vy += particles[i].vy /com.m * particles[i+index].m; 
                    com_shift.vz += particles[i].vz /com.m * particles[i+index].m; 
                    
                    com_shift.x  -= particles[i].x  * particles[i+index_1st_order_a].m/com.m/com.m*dmb; 
                    com_shift.y  -= particles[i].y  * particles[i+index_1st_order_a].m/com.m/com.m*dmb; 
                    com_shift.z  -= particles[i].z  * particles[i+index_1st_order_a].m/com.m/com.m*dmb; 
                    com_shift.vx -= particles[i].vx * particles[i+index_1st_order_a].m/com.m/com.m*dmb; 
                    com_shift.vy -= particles[i].vy * particles[i+index_1st_order_a].m/com.m/com.m*dmb; 
                    com_shift.vz -= particles[i].vz * particles[i+index_1st_order_a].m/com.m/com.m*dmb; 
                    
                    com_shift.x  -= particles[i+index_1st_order_b].x  * particles[i].m/com.m/com.m*dma; 
                    com_shift.y  -= particles[i+index_1st_order_b].y  * particles[i].m/com.m/com.m*dma; 
                    com_shift.z  -= particles[i+index_1st_order_b].z  * particles[i].m/com.m/com.m*dma; 
                    com_shift.vx -= particles[i+index_1st_order_b].vx * particles[i].m/com.m/com.m*dma; 
                    com_shift.vy -= particles[i+index_1st_order_b].vy * particles[i].m/com.m/com.m*dma; 
                    com_shift.vz -= particles[i+index_1st_order_b].vz * particles[i].m/com.m/com.m*dma; 
                    
                    com_shift.x  -= particles[i].x  * particles[i+index_1st_order_b].m/com.m/com.m*dma; 
                    com_shift.y  -= particles[i].y  * particles[i+index_1st_order_b].m/com.m/com.m*dma; 
                    com_shift.z  -= particles[i].z  * particles[i+index_1st_order_b].m/com.m/com.m*dma; 
                    com_shift.vx -= particles[i].vx * particles[i+index_1st_order_b].m/com.m/com.m*dma; 
                    com_shift.vy -= particles[i].vy * particles[i+index_1st_order_b].m/com.m/com.m*dma; 
                    com_shift.vz -= particles[i].vz * particles[i+index_1st_order_b].m/com.m/com.m*dma; 
                    
                    com_shift.x  += 2.*particles[i].x  * particles[i].m/com.m/com.m/com.m*dma*dmb; 
                    com_shift.y  += 2.*particles[i].y  * particles[i].m/com.m/com.m/com.m*dma*dmb; 
                    com_shift.z  += 2.*particles[i].z  * particles[i].m/com.m/com.m/com.m*dma*dmb; 
                    com_shift.vx += 2.*particles[i].vx * particles[i].m/com.m/com.m/com.m*dma*dmb; 
                    com_shift.vy += 2.*particles[i].vy * particles[i].m/com.m/com.m/com.m*dma*dmb; 
                    com_shift.vz += 2.*particles[i].vz * particles[i].m/com.m/com.m/com.m*dma*dmb; 
                    
                    com_shift.x  -= particles[i].x  * particles[i].m/com.m/com.m*ddm; 
                    com_shift.y  -= particles[i].y  * particles[i].m/com.m/com.m*ddm; 
                    com_shift.z  -= particles[i].z  * particles[i].m/com.m/com.m*ddm; 
                    com_shift.vx -= particles[i].vx * particles[i].m/com.m/com.m*ddm; 
                    com_shift.vy -= particles[i].vy * particles[i].m/com.m/com.m*ddm; 
                    com_shift.vz -= particles[i].vz * particles[i].m/com.m/com.m*ddm; 
                    
                    
                    
                    
                    
                }
                for (int i=0;i<N_real;i++){
                    particles[i+index].x -= com_shift.x; 
                    particles[i+index].y -= com_shift.y; 
                    particles[i+index].z -= com_shift.z; 
                    particles[i+index].vx -= com_shift.vx; 
                    particles[i+index].vy -= com_shift.vy; 
                    particles[i+index].vz -= com_shift.vz; 
                }
            }
        }
    }
    // Then do first order
    for (int v=0;v<r->var_config_N;v++){
        int index = r->var_config[v].index;
        if (r->var_config[v].testparticle>=0){
            // Test particles do not affect the COM
        }else{
            if (r->var_config[v].order==1){
                struct reb_particle com_shift = {0};
                double dm = 0.;
                for (int i=0;i<N_real;i++){
                    dm += particles[i+index].m;
                }
                for (int i=0;i<N_real;i++){
                    com_shift.x  += particles[i].m/com.m * particles[i+index].x ; 
                    com_shift.y  += particles[i].m/com.m * particles[i+index].y ; 
                    com_shift.z  += particles[i].m/com.m * particles[i+index].z ; 
                    com_shift.vx += particles[i].m/com.m * particles[i+index].vx; 
                    com_shift.vy += particles[i].m/com.m * particles[i+index].vy; 
                    com_shift.vz += particles[i].m/com.m * particles[i+index].vz; 
                    
                    com_shift.x  += particles[i].x /com.m * particles[i+index].m; 
                    com_shift.y  += particles[i].y /com.m * particles[i+index].m; 
                    com_shift.z  += particles[i].z /com.m * particles[i+index].m; 
                    com_shift.vx += particles[i].vx/com.m * particles[i+index].m; 
                    com_shift.vy += particles[i].vy/com.m * particles[i+index].m; 
                    com_shift.vz += particles[i].vz/com.m * particles[i+index].m; 
                    
                    com_shift.x  -= particles[i].x /(com.m*com.m) * particles[i].m*dm; 
                    com_shift.y  -= particles[i].y /(com.m*com.m) * particles[i].m*dm; 
                    com_shift.z  -= particles[i].z /(com.m*com.m) * particles[i].m*dm; 
                    com_shift.vx -= particles[i].vx/(com.m*com.m) * particles[i].m*dm; 
                    com_shift.vy -= particles[i].vy/(com.m*com.m) * particles[i].m*dm; 
                    com_shift.vz -= particles[i].vz/(com.m*com.m) * particles[i].m*dm; 
                }
                for (int i=0;i<N_real;i++){
                    particles[i+index].x -= com_shift.x; 
                    particles[i+index].y -= com_shift.y; 
                    particles[i+index].z -= com_shift.z; 
                    particles[i+index].vx -= com_shift.vx; 
                    particles[i+index].vy -= com_shift.vy; 
                    particles[i+index].vz -= com_shift.vz; 
                }
            }
        }
    }
	
    // Finally do normal particles
    for (int i=0;i<N_real;i++){
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
    const int N_real = r->N - r->N_var;
	struct reb_particle* restrict const particles = r->particles;
	for (int i=0;i<N_real;i++){
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
		double F = E - e*sin(E) - M;
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
		E = M/fabs(M)*log(2.*fabs(M)/e + 1.8);

		double F = E - e*sinh(E) + M;
		for(int i=0; i<100; i++){
			E = E - F/(1.0 - e*cosh(E));
			F = E - e*sinh(E) + M;
			if(fabs(F) < 1.e-16){
				break;
			}
		}
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

static struct reb_particle reb_particle_nan(void){
    struct reb_particle p;
    p.x = nan("");
    p.y = nan("");
    p.z = nan("");
    p.vx = nan("");
    p.vy = nan("");
    p.vz = nan("");
    p.ax = nan("");
    p.ay = nan("");
    p.az = nan("");
    p.m = nan("");
    p.r = nan("");
    p.lastcollision = nan("");
    p.c = NULL;
    p.id = -1;
    p.ap = NULL;
    p.sim = NULL;

    return p;
}

struct reb_particle reb_tools_orbit_to_particle_err(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f, int* err){
	if(e == 1.){
		*err = 1; 		// Can't initialize a radial orbit with orbital elements.
		return reb_particle_nan();
	}
	if(e < 0.){
		*err = 2; 		// Eccentricity must be greater than or equal to zero.
		return reb_particle_nan();
	}
	if(e > 1.){
		if(a > 0.){
			*err = 3; 	// Bound orbit (a > 0) must have e < 1. 
			return reb_particle_nan();
		}
	}
	else{
		if(a < 0.){
			*err =4; 	// Unbound orbit (a < 0) must have e > 1.
			return reb_particle_nan();
		}
	}
	if(e*cos(f) < -1.){
		*err = 5;		// Unbound orbit can't have f set beyond the range allowed by the asymptotes set by the parabola.
		return reb_particle_nan();
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

struct reb_orbit reb_orbit_nan(void){
    struct reb_orbit o;
    o.d = nan("");
    o.v = nan("");
    o.h = nan("");
    o.P = nan("");
    o.n = nan("");
    o.a = nan("");
    o.e = nan("");
    o.inc = nan("");
    o.Omega = nan("");
    o.omega = nan("");
    o.pomega = nan("");
    o.f = nan("");
    o.M = nan("");
    o.l = nan("");
    o.theta = nan("");
    o.T = nan("");

    return o;
}

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
		return reb_orbit_nan();
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
	o.d = sqrt ( dx*dx + dy*dy + dz*dz );
	
	vsquared = dvx*dvx + dvy*dvy + dvz*dvz;
	o.v = sqrt(vsquared);
	vcircsquared = mu/o.d;	
	o.a = -mu/( vsquared - 2.*vcircsquared );	// semi major axis
	
	hx = (dy*dvz - dz*dvy); 					//angular momentum vector
	hy = (dz*dvx - dx*dvz);
	hz = (dx*dvy - dy*dvx);
	o.h = sqrt ( hx*hx + hy*hy + hz*hz );		// abs value of angular momentum

	vdiffsquared = vsquared - vcircsquared;	
	if(o.d <= TINY){							
		*err = 2;									// particle is on top of primary
		return reb_orbit_nan();
	}
	vr = (dx*dvx + dy*dvy + dz*dvz)/o.d;	
	rvr = o.d*vr;
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

    if(o.e < 1.){
	    ea = acos2(1.-o.d/o.a, o.e, vr);// from definition of eccentric anomaly.  If vr < 0, must be going from apo to peri, so ea = [pi, 2pi] so ea = -acos(cosea)
	    o.M = ea - o.e*sin(ea);			// mean anomaly (Kepler's equation)
    }
    else{
        ea = acosh((1.-o.d/o.a)/o.e);
        if(vr < 0.){                    // Approaching pericenter, so eccentric anomaly < 0.
            ea = -ea;
        }
        o.M = o.e*sinh(ea) - ea;          // Hyperbolic Kepler's equation
    }

	// in the near-planar case, the true longitude is always well defined for the position, and pomega for the pericenter if e!= 0
	// we therefore calculate those and calculate the remaining angles from them
	if(o.inc < MIN_INC || o.inc > M_PI - MIN_INC){	// nearly planar.  Use longitudes rather than angles referenced to node for numerical stability.
		o.pomega = acos2(ex, o.e, ey);		// cos pomega is dot product of x and e unit vectors.  Will = 0 if e=0.
		o.theta = acos2(dx, o.d, dy);		// cos theta is dot product of x and r vectors (true longitude).  Will = 0 if e = 0.
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
		double wpf = acos2(nx*dx + ny*dy, n*o.d, dz);	// omega plus f.  Both angles measured in orbital plane, and always well defined for i!=0.
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
    
    if (p.sim == NULL){                         // if particle isn't in simulation yet, can't get time.  You can still manually apply the equation below using o.M and o.n
        o.T = nan("");
    }
    else{
        o.T = p.sim->t - o.M/fabs(o.n);         // time of pericenter passage (M = n(t-T).  Works for hyperbolic with fabs and n defined as above).
    }

	return o;
}


struct reb_orbit reb_tools_particle_to_orbit(double G, struct reb_particle p, struct reb_particle primary){
	int err;
	return reb_tools_particle_to_orbit_err(G, p, primary, &err);
}

/***********************************
 * Variational Equations and Megno */


int reb_add_var_1st_order(struct reb_simulation* const r, int testparticle){
    r->var_config_N++;
    r->var_config = realloc(r->var_config,sizeof(struct reb_variational_configuration)*r->var_config_N);
    r->var_config[r->var_config_N-1].sim = r;
    r->var_config[r->var_config_N-1].order = 1;
    int index = r->N;
    r->var_config[r->var_config_N-1].index = index;
    r->var_config[r->var_config_N-1].testparticle = testparticle;
    struct reb_particle p0 = {0};
    if (testparticle>=0){
        reb_add(r,p0);
        r->N_var++;
    }else{
        int N_real = r->N - r->N_var;
        for (int i=0;i<N_real;i++){
            reb_add(r,p0);
        }
        r->N_var += N_real;
    }
    return index;
}


int reb_add_var_2nd_order(struct reb_simulation* const r, int testparticle, int index_1st_order_a, int index_1st_order_b){
    r->var_config_N++;
    r->var_config = realloc(r->var_config,sizeof(struct reb_variational_configuration)*r->var_config_N);
    r->var_config[r->var_config_N-1].sim = r;
    r->var_config[r->var_config_N-1].order = 2;
    int index = r->N;
    r->var_config[r->var_config_N-1].index = index;
    r->var_config[r->var_config_N-1].testparticle = testparticle;
    r->var_config[r->var_config_N-1].index_1st_order_a = index_1st_order_a;
    r->var_config[r->var_config_N-1].index_1st_order_b = index_1st_order_b;
    struct reb_particle p0 = {0};
    if (testparticle>=0){
        reb_add(r,p0);
        r->N_var++;
    }else{
        int N_real = r->N - r->N_var;
        for (int i=0;i<N_real;i++){
            reb_add(r,p0);
        }
        r->N_var += N_real;
    }
    return index;
}

#ifndef LIBREBOUNDX
void reb_tools_megno_init(struct reb_simulation* const r){
	r->megno_Ys = 0.;
	r->megno_Yss = 0.;
	r->megno_cov_Yt = 0.;
	r->megno_var_t = 0.;
	r->megno_n = 0;
	r->megno_mean_Y = 0;
	r->megno_mean_t = 0;
    int i = reb_add_var_1st_order(r,-1);
	r->calculate_megno = i;
    const int imax = i + (r->N-r->N_var);
    struct reb_particle* const particles = r->particles;
    for (;i<imax;i++){ 
        particles[i].m  = 0.;
		particles[i].x  = reb_random_normal(1.);
		particles[i].y  = reb_random_normal(1.);
		particles[i].z  = reb_random_normal(1.);
		particles[i].vx = reb_random_normal(1.);
		particles[i].vy = reb_random_normal(1.);
		particles[i].vz = reb_random_normal(1.);
		double deltad = 1./sqrt(
                particles[i].x*particles[i].x 
                + particles[i].y*particles[i].y 
                + particles[i].z*particles[i].z 
                + particles[i].vx*particles[i].vx 
                + particles[i].vy*particles[i].vy 
                + particles[i].vz*particles[i].vz); // rescale
		particles[i].x *= deltad;
		particles[i].y *= deltad;
		particles[i].z *= deltad;
		particles[i].vx *= deltad;
		particles[i].vy *= deltad;
		particles[i].vz *= deltad;
    }
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
    double deltad = 0;
    double delta2 = 0;
    int i = r->calculate_megno;
    const int imax = i + (r->N-r->N_var);
    for (;i<imax;i++){
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

/**************************************
 * New derivatives */

void reb_solve_kepler_pal(double h, double k, double lambda, double* p, double* q){
    double e2 = h*h + k*k;
    if (e2<0.3*0.3){ // low e case
        double pn = 0;
        double qn = 0;

        int n=0;
        double f=0.;
        do{
            double f0 = qn*cos(pn)+pn*sin(pn)-(k*cos(lambda)+h*sin(lambda));
            double f1 = -qn*sin(pn)+pn*cos(pn)-(k*sin(lambda)-h*cos(lambda));

            double fac = 1./(qn-1.);
            double fd00 = fac*(qn*cos(pn)-cos(pn)+pn*sin(pn));
            double fd01 = fac*(pn*cos(pn)-qn*sin(pn)+sin(pn));
            double fd10 = fac*(-sin(pn));
            double fd11 = fac*(-cos(pn));

            qn -= fd00*f0+fd10*f1;
            pn -= fd01*f0+fd11*f1;
            f = sqrt(f0*f0+f1*f1);
        }while(n++<50 && f>1e-15);
        *p = pn;
        *q = qn;
    }else{  // high e case
        double pomega = atan2(h,k);
        double M = lambda-pomega;
        double e = sqrt(e2);
        double E = reb_M_to_E(e, M);
        *p = e*sin(E); 
        *q = e*cos(E); 
    }
}

void reb_particle_to_pal(double G, struct reb_particle p, struct reb_particle primary, double *a, double* lambda, double* k, double* h, double* ix, double* iy){
    double x = p.x - primary.x;
    double y = p.y - primary.y;
    double z = p.z - primary.z;
    double vx = p.vx - primary.vx;
    double vy = p.vy - primary.vy;
    double vz = p.vz - primary.vz;
    double mu = G*(p.m+primary.m);
    double r2 = x*x + y*y + z*z;
    double r = sqrt(r2);
    double cx = y*vz - z*vy;
    double cy = z*vx - x*vz;
    double cz = x*vy - y*vx;
    double c2 = cx*cx + cy*cy + cz*cz;
    double c = sqrt(c2);
    double chat = x*vx + y*vy + z*vz;

    double fac = sqrt(2./(1.+cz/c))/c;
    *ix = -fac * cy;
    *iy = fac * cx;
    *k = c/mu*(vy-vz/(c+cz)*cy)-1./r*(x-z/(c+cz)*cx);
    *h = c/mu*(-vx+vz/(c+cz)*cx)-1./r*(y-z/(c+cz)*cy);
    double e2 = (*k)*(*k)+(*h)*(*h);
    *a = c2/(mu*(1.-e2));
    double l = 1.-sqrt(1.-e2);
    *lambda = atan2(-r*vx+r*vz*cx/(c+cz)-(*k)*chat/(2.-l), r*vy-r*vz*cy/(c+cz)+(*h)*chat/(2.-l))-chat/c*(1.-l);
}

struct reb_particle reb_pal_to_particle(double G, struct reb_particle primary, double m, double a, double lambda, double k, double h, double ix, double iy){
    struct reb_particle np = {0.};
    np.m = m;

    double p=0.,q=0.;
    reb_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double xi = a*clp + p/(2.-l)*h -k;
    double eta = a*slp - p/(2.-l)*k -h;

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double W = eta*ix-xi*iy;

    np.x = primary.x + xi+0.5*iy*W;
    np.y = primary.y + eta-0.5*ix*W;
    np.z = primary.z + 0.5*iz*W;

    double an = sqrt(G*(m+primary.m)/a);
    double dxi  = an/(1.-q)*(-slp+q/(2.-l)*h);
    double deta = an/(1.-q)*(+clp-q/(2.-l)*k);
    double dW = deta*ix-dxi*iy;

    np.vx = primary.vx + dxi+0.5*iy*dW;
    np.vy = primary.vy + deta-0.5*ix*dW;
    np.vz = primary.vz + 0.5*iz*dW;


    return np;
}


struct reb_particle reb_vary_pal_lambda(double G, struct reb_particle po, struct reb_particle primary){
    double a, lambda, k, h, ix, iy;
    reb_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};

    double p=0.,q=0.;
    reb_solve_kepler_pal(h, k, lambda, &p, &q);
    double dq_dlambda = -p/(1.-q);
    double dp_dlambda = q/(1.-q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dlambda = -1./(1.-q)*slp;
    double dslp_dlambda = 1./(1.-q)*clp;
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dxi_dlambda = a*dclp_dlambda + dp_dlambda/(2.-l)*h;
    double deta_dlambda = a*dslp_dlambda - dp_dlambda/(2.-l)*k;

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dlambda = deta_dlambda*ix-dxi_dlambda*iy;

    np.x = dxi_dlambda+0.5*iy*dW_dlambda;
    np.y = deta_dlambda-0.5*ix*dW_dlambda;
    np.z = 0.5*iz*dW_dlambda;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dlambda  = an/((1.-q)*(1.-q))*dq_dlambda*(-slp+q/(2.-l)*h)
    + an/(1.-q)*(-dslp_dlambda+dq_dlambda/(2.-l)*h);
    double ddeta_dlambda = an/((1.-q)*(1.-q))*dq_dlambda*(+clp-q/(2.-l)*k)
    + an/(1.-q)*(dclp_dlambda-dq_dlambda/(2.-l)*k);
    double ddW_dlambda = ddeta_dlambda*ix-ddxi_dlambda*iy;
    np.vx = ddxi_dlambda+0.5*iy*ddW_dlambda;
    np.vy = ddeta_dlambda-0.5*ix*ddW_dlambda;
    np.vz = 0.5*iz*ddW_dlambda;

    return np;
}

struct reb_particle reb_vary_pal_h(double G, struct reb_particle po, struct reb_particle primary){
    double a, lambda, k, h, ix, iy;
    reb_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dh = -1./(1.-q)*(-slp*clp);
    double dslp_dh = -1./(1.-q)*(clp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dh = 1./sqrt(1.-h*h-k*k)*h;
    double dp_dh = 1./(1.-q)*(-clp);
    double dxi_dh = a*dclp_dh + dp_dh/(2.-l)*h + p/(2.-l) + p/((2.-l)*(2.-l))*dl_dh*h;
    double deta_dh = a*dslp_dh - dp_dh/(2.-l)*k - p/((2.-l)*(2.-l))*k*dl_dh -1;

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dh = deta_dh*ix-dxi_dh*iy;

    np.x = dxi_dh+0.5*iy*dW_dh;
    np.y = deta_dh-0.5*ix*dW_dh;
    np.z = 0.5*iz*dW_dh;

    double dq_dh = 1./(1.-q)*(slp-h);

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dh  = dq_dh*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h)
                + an/(1.-q) * (-dslp_dh+dq_dh/(2.-l)*h+dl_dh*q/((2.-l)*(2.-l))*h+q/(2.-l));
    double ddeta_dh = dq_dh*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k)
                + an/(1.-q) * (+dclp_dh-dq_dh/(2.-l)*k-dl_dh*q/((2.-l)*(2.-l))*k);
    double ddW_dh = ddeta_dh*ix-ddxi_dh*iy;

    np.vx = ddxi_dh+0.5*iy*ddW_dh;
    np.vy = ddeta_dh-0.5*ix*ddW_dh;
    np.vz = 0.5*iz*ddW_dh;

    return np;
}

struct reb_particle reb_vary_pal_k(double G, struct reb_particle po, struct reb_particle primary){
    double a, lambda, k, h, ix, iy;
    reb_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dk = -1./(1.-q)*(slp*slp);
    double dslp_dk = -1./(1.-q)*(-slp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dk = 1./sqrt(1.-h*h-k*k)*k;
    double dp_dk = 1./(1.-q)*(slp);
    double dxi_dk = a*dclp_dk + dp_dk/(2.-l)*h + p/((2.-l)*(2.-l))*dl_dk*h -1;
    double deta_dk = a*dslp_dk - dp_dk/(2.-l)*k - p/(2.-l) - p/((2.-l)*(2.-l))*dl_dk*k;

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dk = deta_dk*ix-dxi_dk*iy;

    np.x = dxi_dk+0.5*iy*dW_dk;
    np.y = deta_dk-0.5*ix*dW_dk;
    np.z = 0.5*iz*dW_dk;

    double dq_dk = 1./(1.-q)*(clp-k);

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dk  = dq_dk*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h)
                + an/(1.-q) * (-dslp_dk+dq_dk/(2.-l)*h+dl_dk*q/((2.-l)*(2.-l))*h);
    double ddeta_dk = dq_dk*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k)
                + an/(1.-q) * (+dclp_dk-dq_dk/(2.-l)*k-dl_dk*q/((2.-l)*(2.-l))*k-q/(2.-l));
    double ddW_dk = ddeta_dk*ix-ddxi_dk*iy;

    np.vx = ddxi_dk+0.5*iy*ddW_dk;
    np.vy = ddeta_dk-0.5*ix*ddW_dk;
    np.vz = 0.5*iz*ddW_dk;

    return np;
}

struct reb_particle reb_vary_pal_kk(double G, struct reb_particle po, struct reb_particle primary){
    double a, lambda, k, h, ix, iy;
    reb_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dk = -1./(1.-q)*(slp*slp);
    double dslp_dk = -1./(1.-q)*(-slp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dk = 1./sqrt(1.-h*h-k*k)*k;
    double dl_dkk = 1./sqrt(1.-h*h-k*k) + (k*k)/(sqrt(1.-h*h-k*k)*sqrt(1.-h*h-k*k)*sqrt(1.-h*h-k*k));
    double dp_dk = 1./(1.-q)*(slp);
    double dq_dk = 1./(1.-q)*(clp-k);
    double dp_dkk = dq_dk/((1.-q)*(1.-q))*(slp) + 1./(1.-q)*(dslp_dk);
    double dq_dkk = dq_dk/((1.-q)*(1.-q))*(clp-k) + 1./(1.-q)*(dclp_dk -1.);
    double dclp_dkk = -dq_dk/((1.-q)*(1.-q))*(slp*slp) -2./(1.-q)*slp*dslp_dk;
    double dslp_dkk = -dq_dk/((1.-q)*(1.-q))*(-slp*clp) -1./(1.-q)*-slp*dclp_dk -1./(1.-q)*-dslp_dk*clp;

    double dxi_dkk = a*dclp_dkk + dp_dkk/(2.-l)*h + dl_dk*dp_dk/((2.-l)*(2.-l))*h + dp_dk/((2.-l)*(2.-l))*dl_dk*h + 2.*dl_dk*p/((2.-l)*(2.-l)*(2.-l))*dl_dk*h + p/((2.-l)*(2.-l))*dl_dkk*h;
    double deta_dkk = a*dslp_dkk - dp_dkk/(2.-l)*k - dl_dk*dp_dk/((2.-l)*(2.-l))*k - dp_dk/(2.-l) - dp_dk/(2.-l) - dl_dk*p/((2.-l)*(2.-l)) 
                - dp_dk/((2.-l)*(2.-l))*dl_dk*k - 2.*dl_dk*p/((2.-l)*(2.-l)*(2.-l))*dl_dk*k - p/((2.-l)*(2.-l))*dl_dkk*k - p/((2.-l)*(2.-l))*dl_dk;

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dkk = deta_dkk*ix-dxi_dkk*iy;

    np.x = dxi_dkk+0.5*iy*dW_dkk;
    np.y = deta_dkk-0.5*ix*dW_dkk;
    np.z = 0.5*iz*dW_dkk;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dkk  = dq_dkk*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h) + 2.*dq_dk*dq_dk*an/((1.-q)*(1.-q)*(1.-q))*(-slp+q/(2.-l)*h)
    			+ dq_dk*an/((1.-q)*(1.-q))*(-dslp_dk+dq_dk/(2.-l)*h+dl_dk*q/((2.-l)*(2.-l))*h)
                + dq_dk*an/((1.-q)*(1.-q))*(-dslp_dk+dq_dk/(2.-l)*h+dl_dk*q/((2.-l)*(2.-l))*h)
                + an/(1.-q)*(-dslp_dkk + dq_dkk/(2.-l)*h + dl_dk*dq_dk/((2.-l)*(2.-l))*h 
                + dl_dkk*q/((2.-l)*(2.-l))*h + dl_dk*dq_dk/((2.-l)*(2.-l))*h + 2.*dl_dk*dl_dk*q/((2.-l)*(2.-l)*(2.-l))*h );
    double ddeta_dkk = dq_dkk*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k) + 2.*dq_dk*dq_dk*an/((1.-q)*(1.-q)*(1.-q))*(+clp-q/(2.-l)*k) 
    			+ dq_dk*an/((1.-q)*(1.-q))*(+dclp_dk-dq_dk/(2.-l)*k-dl_dk*q/((2.-l)*(2.-l))*k-q/(2.-l))
                + dq_dk*an/((1.-q)*(1.-q))*(+dclp_dk-dq_dk/(2.-l)*k-dl_dk*q/((2.-l)*(2.-l))*k-q/(2.-l))
                + an/(1.-q)*(+dclp_dkk - dq_dkk/(2.-l)*k - dq_dk*dl_dk/((2.-l)*(2.-l))*k - dq_dk/(2.-l) 
                - dl_dkk*q/((2.-l)*(2.-l))*k - dl_dk*dq_dk/((2.-l)*(2.-l))*k - 2.*dl_dk*dl_dk*q/((2.-l)*(2.-l)*(2.-l))*k - dl_dk*q/((2.-l)*(2.-l)) - dq_dk/(2.-l) - dl_dk*q/((2.-l)*(2.-l)) );
    double ddW_dkk = ddeta_dkk*ix-ddxi_dkk*iy;

    np.vx = ddxi_dkk+0.5*iy*ddW_dkk;
    np.vy = ddeta_dkk-0.5*ix*ddW_dkk;
    np.vz = 0.5*iz*ddW_dkk;

    return np;
}

struct reb_particle reb_vary_pal_hh(double G, struct reb_particle po, struct reb_particle primary){
    double a, lambda, k, h, ix, iy;
    reb_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dh = -1./(1.-q)*(-slp*clp);
    double dslp_dh = -1./(1.-q)*(clp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dh = 1./sqrt(1.-h*h-k*k)*h;
    double dl_dhh = 1./sqrt(1.-h*h-k*k) + (h*h)/(sqrt(1.-h*h-k*k)*sqrt(1.-h*h-k*k)*sqrt(1.-h*h-k*k));
    double dp_dh = 1./(1.-q)*(-clp);
    double dq_dh = 1./(1.-q)*(slp-h);
	double dq_dhh = 1./((1.-q)*(1.-q))*dq_dh*(slp-h) + 1./(1.-q)*(dslp_dh-1);
	double dp_dhh = 1./((1.-q)*(1.-q))*dq_dh*(-clp) + 1./(1.-q)*(-dclp_dh);
    double dclp_dhh = -1./((1.-q)*(1.-q))*dq_dh*(-slp*clp) - 1./(1.-q)*(-dslp_dh*clp) - 1./(1.-q)*(-slp*dclp_dh);
    double dslp_dhh = -1./((1.-q)*(1.-q))*dq_dh*(clp*clp) - 2./(1.-q)*(clp*dclp_dh);

    double dxi_dhh = a*dclp_dhh + (dp_dhh/(2.-l)*h + dl_dh*dp_dh/((2.-l)*(2.-l))*h + dp_dh/(2.-l)) + (dp_dh/(2.-l)+ dl_dh*p/((2.-l)*(2.-l)))
    	+ (dp_dh/((2.-l)*(2.-l))*dl_dh*h + 2.*p/((2.-l)*(2.-l)*(2.-l))*dl_dh*dl_dh*h + p/((2.-l)*(2.-l))*dl_dhh*h + p/((2.-l)*(2.-l))*dl_dh);
    double deta_dhh = a*dslp_dhh + (-dp_dhh/(2.-l)*k - dl_dh*dp_dh/((2.-l)*(2.-l))*k) +(- dp_dh/((2.-l)*(2.-l))*k*dl_dh - 2.*p/((2.-l)*(2.-l)*(2.-l))*k*dl_dh*dl_dh- p/((2.-l)*(2.-l))*k*dl_dhh );

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dhh = deta_dhh*ix-dxi_dhh*iy;

    np.x = dxi_dhh+0.5*iy*dW_dhh;
    np.y = deta_dhh-0.5*ix*dW_dhh;
    np.z = 0.5*iz*dW_dhh;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dhh  = dq_dhh*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h) + 2.*dq_dh*dq_dh*an/((1.-q)*(1.-q)*(1.-q))*(-slp+q/(2.-l)*h) 
    			+ dq_dh*an/((1.-q)*(1.-q))*(-dslp_dh+dq_dh/(2.-l)*h+dl_dh*q/((2.-l)*(2.-l))*h+q/(2.-l))
                + dq_dh*an/((1.-q)*(1.-q))*(-dslp_dh+dq_dh/(2.-l)*h+dl_dh*q/((2.-l)*(2.-l))*h+q/(2.-l)) 
                + an/(1.-q)*(-dslp_dhh + (dq_dhh/(2.-l)*h+dl_dh*dq_dh/((2.-l)*(2.-l))*h+dq_dh/(2.-l)) 
                + (dl_dhh*q/((2.-l)*(2.-l))*h+dl_dh*dq_dh/((2.-l)*(2.-l))*h+2.*dl_dh*dl_dh*q/((2.-l)*(2.-l)*(2.-l))*h+dl_dh*q/((2.-l)*(2.-l))) + (dq_dh/(2.-l)+dl_dh*q/((2.-l)*(2.-l)))  );
    double ddeta_dhh = dq_dhh*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k) + 2.*dq_dh*dq_dh*an/((1.-q)*(1.-q)*(1.-q))*(+clp-q/(2.-l)*k) 
    			+ dq_dh*an/((1.-q)*(1.-q))*(+dclp_dh-dq_dh/(2.-l)*k-dl_dh*q/((2.-l)*(2.-l))*k)
    			+ dq_dh*an/((1.-q)*(1.-q))*(+dclp_dh-dq_dh/(2.-l)*k-dl_dh*q/((2.-l)*(2.-l))*k)
                + an/(1.-q)*(+dclp_dhh - dq_dhh/(2.-l)*k - dl_dh*dq_dh/((2.-l)*(2.-l))*k 
                - dl_dhh*q/((2.-l)*(2.-l))*k - dl_dh*dq_dh/((2.-l)*(2.-l))*k - 2.*dl_dh*dl_dh*q/((2.-l)*(2.-l)*(2.-l))*k );

    double ddW_dhh = ddeta_dhh*ix-ddxi_dhh*iy;

    np.vx = ddxi_dhh+0.5*iy*ddW_dhh;
    np.vy = ddeta_dhh-0.5*ix*ddW_dhh;
    np.vz = 0.5*iz*ddW_dhh;

    return np;
}

struct reb_particle reb_vary_pal_lambdalambda(double G, struct reb_particle po, struct reb_particle primary){
    double a, lambda, k, h, ix, iy;
    reb_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};

    double p=0.,q=0.;
    reb_solve_kepler_pal(h, k, lambda, &p, &q);
    double dq_dlambda = -p/(1.-q);
    double dp_dlambda = q/(1.-q);
    double dq_dlambdalambda = -dp_dlambda/(1.-q) - p/((1.-q)*(1.-q))*dq_dlambda ;
    double dp_dlambdalambda = dq_dlambda/(1.-q) + q/((1.-q)*(1.-q))*dq_dlambda ;

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dlambda = -1./(1.-q)*slp;
    double dslp_dlambda = 1./(1.-q)*clp;
    double dclp_dlambdalambda = -1./((1.-q)*(1.-q))*dq_dlambda*slp -1./(1.-q)*dslp_dlambda;
    double dslp_dlambdalambda = 1./((1.-q)*(1.-q))*dq_dlambda*clp + 1./(1.-q)*dclp_dlambda;    
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dxi_dlambdalambda = a*dclp_dlambdalambda + dp_dlambdalambda/(2.-l)*h;
    double deta_dlambdalambda = a*dslp_dlambdalambda - dp_dlambdalambda/(2.-l)*k;

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dlambdalambda = deta_dlambdalambda*ix-dxi_dlambdalambda*iy;

    np.x = dxi_dlambdalambda+0.5*iy*dW_dlambdalambda;
    np.y = deta_dlambdalambda-0.5*ix*dW_dlambdalambda;
    np.z = 0.5*iz*dW_dlambdalambda;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dlambdalambda  = 2.*an/((1.-q)*(1.-q)*(1.-q))*dq_dlambda*dq_dlambda*(-slp+q/(2.-l)*h) 
    			+ an/((1.-q)*(1.-q))*dq_dlambdalambda*(-slp+q/(2.-l)*h) + an/((1.-q)*(1.-q))*dq_dlambda*(-dslp_dlambda+dq_dlambda/(2.-l)*h)
    			+ an/((1.-q)*(1.-q))*dq_dlambda*(-dslp_dlambda+dq_dlambda/(2.-l)*h) + an/(1.-q)*(-dslp_dlambdalambda+dq_dlambdalambda/(2.-l)*h);
    double ddeta_dlambdalambda = 2.*an/((1.-q)*(1.-q)*(1.-q))*dq_dlambda*dq_dlambda*(+clp-q/(2.-l)*k) 
    			+ an/((1.-q)*(1.-q))*dq_dlambdalambda*(+clp-q/(2.-l)*k) + an/((1.-q)*(1.-q))*dq_dlambda*(dclp_dlambda-dq_dlambda/(2.-l)*k)
    			+ an/((1.-q)*(1.-q))*dq_dlambda*(dclp_dlambda-dq_dlambda/(2.-l)*k) + an/(1.-q)*(dclp_dlambdalambda-dq_dlambdalambda/(2.-l)*k);

    double ddW_dlambdalambda = ddeta_dlambdalambda*ix-ddxi_dlambdalambda*iy;
    np.vx = ddxi_dlambdalambda+0.5*iy*ddW_dlambdalambda;
    np.vy = ddeta_dlambdalambda-0.5*ix*ddW_dlambdalambda;
    np.vz = 0.5*iz*ddW_dlambdalambda;

    return np;
}

struct reb_particle reb_vary_pal_klambda(double G, struct reb_particle po, struct reb_particle primary){
    double a, lambda, k, h, ix, iy;
    reb_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dk = -1./(1.-q)*(slp*slp);
    double dslp_dk = -1./(1.-q)*(-slp*clp);
    double dclp_dlambda = -1./(1.-q)*slp;
    double dslp_dlambda = 1./(1.-q)*clp;
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dk = 1./sqrt(1.-h*h-k*k)*k;
    double dp_dk = 1./(1.-q)*(slp);
    double dq_dk = 1./(1.-q)*(clp-k);
    double dq_dlambda = -p/(1.-q);
    double dp_dlambda = q/(1.-q);
    double dq_dklambda = -dp_dk/(1.-q) -p/((1.-q)*(1.-q))*dq_dk;
    double dp_dklambda = dq_dk/(1.-q) + q/((1.-q)*(1.-q))*dq_dk;
    double dclp_dklambda = -1./(1.-q)*dslp_dk -1./((1.-q)*(1.-q))*dq_dk*slp;
    double dslp_dklambda = 1./(1.-q)*dclp_dk + 1./((1.-q)*(1.-q))*dq_dk*clp;


    double dxi_dklambda = a*dclp_dklambda + dp_dklambda/(2.-l)*h + dp_dlambda/((2.-l)*(2.-l))*dl_dk*h;
    double deta_dklambda = a*dslp_dklambda - dp_dklambda/(2.-l)*k - dp_dlambda/(2.-l) - dp_dlambda/((2.-l)*(2.-l))*dl_dk*k;

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dklambda = deta_dklambda*ix-dxi_dklambda*iy;

    np.x = dxi_dklambda+0.5*iy*dW_dklambda;
    np.y = deta_dklambda-0.5*ix*dW_dklambda;
    np.z = 0.5*iz*dW_dklambda;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dklambda  = dq_dklambda*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h) + 2.*dq_dk*dq_dlambda*an/((1.-q)*(1.-q)*(1.-q))*(-slp+q/(2.-l)*h) 
    			+ dq_dk*an/((1.-q)*(1.-q))*(-dslp_dlambda+dq_dlambda/(2.-l)*h)
                + dq_dlambda*an/((1.-q)*(1.-q))*(-dslp_dk+dq_dk/(2.-l)*h+dl_dk*q/((2.-l)*(2.-l))*h)
                + an/(1.-q)*(-dslp_dklambda+dq_dklambda/(2.-l)*h+dl_dk*dq_dlambda/((2.-l)*(2.-l))*h);
    double ddeta_dklambda = dq_dklambda*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k) + 2.*dq_dk*dq_dlambda*an/((1.-q)*(1.-q)*(1.-q))*(+clp-q/(2.-l)*k) 
    			+ dq_dk*an/((1.-q)*(1.-q))*(+dclp_dlambda-dq_dlambda/(2.-l)*k)
                + dq_dlambda*an/((1.-q)*(1.-q))*(+dclp_dk-dq_dk/(2.-l)*k-dl_dk*q/((2.-l)*(2.-l))*k-q/(2.-l)) 
                + an/(1.-q)*(+dclp_dklambda-dq_dklambda/(2.-l)*k-dl_dk*dq_dlambda/((2.-l)*(2.-l))*k-dq_dlambda/(2.-l));
    double ddW_dklambda = ddeta_dklambda*ix-ddxi_dklambda*iy;

    np.vx = ddxi_dklambda+0.5*iy*ddW_dklambda;
    np.vy = ddeta_dklambda-0.5*ix*ddW_dklambda;
    np.vz = 0.5*iz*ddW_dklambda;

    return np;
}

struct reb_particle reb_vary_pal_hlambda(double G, struct reb_particle po, struct reb_particle primary){
    double a, lambda, k, h, ix, iy;
    reb_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dh = -1./(1.-q)*(-slp*clp);
    double dslp_dh = -1./(1.-q)*(clp*clp);
    double dclp_dlambda = -1./(1.-q)*slp;
    double dslp_dlambda = 1./(1.-q)*clp;
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dh = 1./sqrt(1.-h*h-k*k)*h;
    double dp_dh = 1./(1.-q)*(-clp);
    double dq_dh = 1./(1.-q)*(slp-h);
    double dq_dlambda = -p/(1.-q);
    double dp_dlambda = q/(1.-q);
    double dq_dhlambda = -dp_dh/(1.-q) -p/((1.-q)*(1.-q))*dq_dh;
    double dp_dhlambda = dq_dh/(1.-q) +q/((1.-q)*(1.-q))*dq_dh;
    double dclp_dhlambda = -1./((1.-q)*(1.-q))*(-slp*clp)*dq_dlambda -1./(1.-q)*(-dslp_dlambda*clp) -1./(1.-q)*(-slp*dclp_dlambda);
    double dslp_dhlambda = -1./((1.-q)*(1.-q))*(clp*clp)*dq_dlambda -2./(1.-q)*(clp*dclp_dlambda);

    double dxi_dhlambda = a*dclp_dhlambda + dp_dhlambda/(2.-l)*h + dp_dlambda/(2.-l) + dp_dlambda/((2.-l)*(2.-l))*dl_dh*h;
    double deta_dhlambda = a*dslp_dhlambda - dp_dhlambda/(2.-l)*k - dp_dlambda/((2.-l)*(2.-l))*k*dl_dh;

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dhlambda = deta_dhlambda*ix-dxi_dhlambda*iy;

    np.x = dxi_dhlambda+0.5*iy*dW_dhlambda;
    np.y = deta_dhlambda-0.5*ix*dW_dhlambda;
    np.z = 0.5*iz*dW_dhlambda;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dhlambda  = dq_dhlambda*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h) + 2.*dq_dlambda*dq_dh*an/((1.-q)*(1.-q)*(1.-q))*(-slp+q/(2.-l)*h) 
    			+ dq_dh*an/((1.-q)*(1.-q))*(-dslp_dlambda+dq_dlambda/(2.-l)*h)
                + dq_dlambda*an/((1.-q)*(1.-q))*(-dslp_dh+dq_dh/(2.-l)*h+dl_dh*q/((2.-l)*(2.-l))*h+q/(2.-l)) 
                + an/(1.-q)*(-dslp_dhlambda+dq_dhlambda/(2.-l)*h+dl_dh*dq_dlambda/((2.-l)*(2.-l))*h+dq_dlambda/(2.-l));
    double ddeta_dhlambda = dq_dhlambda*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k) + 2.*dq_dh*dq_dlambda*an/((1.-q)*(1.-q)*(1.-q))*(+clp-q/(2.-l)*k) 
    			+ dq_dh*an/((1.-q)*(1.-q))*(+dclp_dlambda-dq_dlambda/(2.-l)*k)
                + dq_dlambda*an/((1.-q)*(1.-q))*(+dclp_dh-dq_dh/(2.-l)*k-dl_dh*q/((2.-l)*(2.-l))*k) 
                + an/(1.-q)*(+dclp_dhlambda-dq_dhlambda/(2.-l)*k-dl_dh*dq_dlambda/((2.-l)*(2.-l))*k);
    double ddW_dhlambda = ddeta_dhlambda*ix-ddxi_dhlambda*iy;

    np.vx = ddxi_dhlambda+0.5*iy*ddW_dhlambda;
    np.vy = ddeta_dhlambda-0.5*ix*ddW_dhlambda;
    np.vz = 0.5*iz*ddW_dhlambda;

    return np;
}

struct reb_particle reb_vary_pal_kh(double G, struct reb_particle po, struct reb_particle primary){
    double a, lambda, k, h, ix, iy;
    reb_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dk = -1./(1.-q)*(slp*slp);
    double dslp_dk = -1./(1.-q)*(-slp*clp);
    double dclp_dh = -1./(1.-q)*(-slp*clp);
    double dslp_dh = -1./(1.-q)*(clp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dh = 1./sqrt(1.-h*h-k*k)*h;
    double dp_dh = 1./(1.-q)*(-clp);
    double dq_dh = 1./(1.-q)*(slp-h);
    double dl_dk = 1./sqrt(1.-h*h-k*k)*k;
    double dp_dk = 1./(1.-q)*(slp);
    double dq_dk = 1./(1.-q)*(clp-k);
    double dl_dkh = k*h/(sqrt(1.-h*h-k*k)*sqrt(1.-h*h-k*k)*sqrt(1.-h*h-k*k));
    double dp_dkh = 1./((1.-q)*(1.-q))*dq_dh*(slp) + 1./(1.-q)*(dslp_dh);
    double dq_dkh = 1./((1.-q)*(1.-q))*dq_dh*(clp-k) + 1./(1.-q)*(dclp_dh);
    double dclp_dkh = -1./((1.-q)*(1.-q))*dq_dh*(slp*slp) -2./(1.-q)*(slp*dslp_dh);
    double dslp_dkh = -1./((1.-q)*(1.-q))*dq_dh*(-slp*clp) -1./(1.-q)*(-dslp_dh*clp) -1./(1.-q)*(-slp*dclp_dh);

    double dxi_dkh = a*dclp_dkh + dp_dkh/(2.-l)*h + dl_dh*dp_dk/((2.-l)*(2.-l))*h + dp_dk/(2.-l) 
    			+ dp_dh/((2.-l)*(2.-l))*dl_dk*h + 2.*p/((2.-l)*(2.-l)*(2.-l))*dl_dk*dl_dh*h + p/((2.-l)*(2.-l))*dl_dkh*h + p/((2.-l)*(2.-l))*dl_dk;
    double deta_dkh = a*dslp_dkh - dp_dkh/(2.-l)*k - dl_dh*dp_dk/((2.-l)*(2.-l))*k - dp_dh/(2.-l)- dl_dh*p/((2.-l)*(2.-l)) 
                - dp_dh/((2.-l)*(2.-l))*dl_dk*k - p/((2.-l)*(2.-l))*dl_dkh*k - 2.*p/((2.-l)*(2.-l)*(2.-l))*dl_dk*dl_dh*k;

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dkh = deta_dkh*ix-dxi_dkh*iy;

    np.x = dxi_dkh+0.5*iy*dW_dkh;
    np.y = deta_dkh-0.5*ix*dW_dkh;
    np.z = 0.5*iz*dW_dkh;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dkh = dq_dkh*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h) + 2.*dq_dh*dq_dk*an/((1.-q)*(1.-q)*(1.-q))*(-slp+q/(2.-l)*h) 
    			+ dq_dk*an/((1.-q)*(1.-q))*(-dslp_dh+dq_dh/(2.-l)*h + dl_dh*q/((2.-l)*(2.-l))*h + q/(2.-l))
                + dq_dh*an/((1.-q)*(1.-q))*(-dslp_dk+dq_dk/(2.-l)*h+dl_dk*q/((2.-l)*(2.-l))*h) 
                + an/(1.-q)*(-dslp_dkh+(dq_dkh/(2.-l)*h+dl_dh*dq_dk/((2.-l)*(2.-l))*h+dq_dk/(2.-l)) 
                + dl_dkh*q/((2.-l)*(2.-l))*h + dl_dk*dq_dh/((2.-l)*(2.-l))*h + 2.*dl_dh*dl_dk*q/((2.-l)*(2.-l)*(2.-l))*h + dl_dk*q/((2.-l)*(2.-l)));
    double ddeta_dkh = dq_dkh*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k) + 2.*dq_dh*dq_dk*an/((1.-q)*(1.-q)*(1.-q))*(+clp-q/(2.-l)*k) 
    			+ dq_dk*an/((1.-q)*(1.-q))*(+dclp_dh-dq_dh/(2.-l)*k-dl_dh*q/((2.-l)*(2.-l))*k)
                + dq_dh*an/((1.-q)*(1.-q))*(+dclp_dk-dq_dk/(2.-l)*k-dl_dk*q/((2.-l)*(2.-l))*k-q/(2.-l))
                + an/(1.-q)*(+dclp_dkh-dq_dkh/(2.-l)*k-dl_dh*dq_dk/((2.-l)*(2.-l))*k 
                -dl_dkh*q/((2.-l)*(2.-l))*k -dl_dk*dq_dh/((2.-l)*(2.-l))*k -2.*dl_dk*dl_dh*q/((2.-l)*(2.-l)*(2.-l))*k -dq_dh/(2.-l)-dl_dh*q/((2.-l)*(2.-l)) );
    double ddW_dkh = ddeta_dkh*ix-ddxi_dkh*iy;

    np.vx = ddxi_dkh+0.5*iy*ddW_dkh;
    np.vy = ddeta_dkh-0.5*ix*ddW_dkh;
    np.vz = 0.5*iz*ddW_dkh;

    return np;
}

struct reb_particle reb_vary_pal_a(double G, struct reb_particle po, struct reb_particle primary){
    double a, lambda, k, h, ix, iy;
    reb_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};

    double p=0.,q=0.;
    reb_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dxi_da = clp;
    double deta_da = slp;

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_da = deta_da*ix-dxi_da*iy;

    np.x = dxi_da+0.5*iy*dW_da;
    np.y = deta_da-0.5*ix*dW_da;
    np.z = 0.5*iz*dW_da;

    double dan_da = sqrt(G*(po.m+primary.m))* -sqrt(1./(a*a*a))/2.;
    double ddxi_da  = dan_da/(1.-q)*q*(-slp+q/(2.-l)*h);
    double ddeta_da = dan_da/(1.-q)*q*(+clp-q/(2.-l)*k);

    double ddW_da = ddeta_da*ix-ddxi_da*iy;
    np.vx = ddxi_da+0.5*iy*ddW_da;
    np.vy = ddeta_da-0.5*ix*ddW_da;
    np.vz = 0.5*iz*ddW_da;

    return np;
}


/**struct reb_particle reb_pal_to_particle(double G, struct reb_particle primary, double m, double a, double lambda, double k, double h, double ix, double iy){
    struct reb_particle np = {0.};
    np.m = m;

    double p=0.,q=0.;
    reb_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double xi = a*clp + p/(2.-l)*h -k;
    double eta = a*slp - p/(2.-l)*k -h;

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double W = eta*ix-xi*iy;

    np.x = primary.x + xi+0.5*iy*W;
    np.y = primary.y + eta-0.5*ix*W;
    np.z = primary.z + 0.5*iz*W;

    double an = sqrt(G*(m+primary.m)/a);
    double dxi  = an/(1.-q)*(-slp+q/(2.-l)*h);
    double deta = an/(1.-q)*(+clp-q/(2.-l)*k);
    double dW = deta*ix-dxi*iy;

    np.vx = primary.vx + dxi+0.5*iy*dW;
    np.vy = primary.vy + deta-0.5*ix*dW;
    np.vz = primary.vz + 0.5*iz*dW;


    return np;
}*/

/**************************************
 * Functionis for derivates of orbits  */

struct reb_particle reb_tools_orbit_to_particle_da(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){

	struct reb_particle p = {0};
	double dr = (1.-e*e)/(1. + e*cos(f));
	double dv0 = -0.5/sqrt(a*a*a)*sqrt(G*(m+primary.m)/(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.x = dr*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
	p.y = dr*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
	p.z = dr*(so*cf+co*sf)*si;

	p.vx = dv0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
	p.vy = dv0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
	p.vz = dv0*((e+cf)*co*si - sf*si*so);
	
	return p;
}

struct reb_particle reb_tools_orbit_to_particle_dda(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){

	struct reb_particle p = {0};
	double ddv0 = 0.75/(a*a*sqrt(a))*sqrt(G*(m+primary.m)/(1.-e*e));

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.vx = ddv0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
	p.vy = ddv0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
	p.vz = ddv0*((e+cf)*co*si - sf*si*so);

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_de(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
    double cosf = cos(f);
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 
	double dr = -a*(cosf*e*e+cosf+2.*e)/((cosf*e+1.)*(cosf*e+1.));
	double dv0 = sqrt(G*(m+primary.m)/a)*e/((1.-e*e)*sqrt(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.x = dr*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
	p.y = dr*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
	p.z = dr*(so*cf+co*sf)*si;

	p.vx = dv0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
	p.vy = dv0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
	p.vz = dv0*((e+cf)*co*si - sf*si*so);
	
    p.vx += v0*(-ci*co*sO - cO*so);
	p.vy += v0*(ci*co*cO - sO*so);
	p.vz += v0*(co*si);

	return p;
}



struct reb_particle reb_tools_orbit_to_particle_dde(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
    double cosf = cos(f);
	double ddr = a*2.*(cosf*cosf-1.)/((cosf*e+1.)*(cosf*e+1.)*(cosf*e+1.));
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 
	double dv0 = e*v0/(1.-e*e); 
	double ddv0 = v0/((e*e-1.)*(e*e-1.)) * (2.*e*e+1.);

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.x = ddr*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
	p.y = ddr*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
	p.z = ddr*(so*cf+co*sf)*si;

	p.vx = ddv0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
	p.vy = ddv0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
	p.vz = ddv0*((e+cf)*co*si - sf*si*so);
	
    p.vx += 2.*dv0*(-ci*co*sO - cO*so);
	p.vy += 2.*dv0*(ci*co*cO - sO*so);
	p.vz += 2.*dv0*(co*si);

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_di(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double r = a*(1.-e*e)/(1. + e*cos(f));
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double dci = -sin(inc);
	double dsi = cos(inc);
	
	p.x = r*(- sO*(so*cf+co*sf)*dci);
	p.y = r*(+ cO*(so*cf+co*sf)*dci);
	p.z = r*(so*cf+co*sf)*dsi;

	p.vx = v0*((e+cf)*(-dci*co*sO) - sf*(- dci*so*sO));
	p.vy = v0*((e+cf)*(dci*co*cO)  - sf*(dci*so*cO));
	p.vz = v0*((e+cf)*co*dsi - sf*dsi*so);
	

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_ddi(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double r = a*(1.-e*e)/(1. + e*cos(f));
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ddci = -cos(inc);
	double ddsi = -sin(inc);
	
	p.x = r*(- sO*(so*cf+co*sf)*ddci);
	p.y = r*(+ cO*(so*cf+co*sf)*ddci);
	p.z = r*(so*cf+co*sf)*ddsi;

	p.vx = v0*((e+cf)*(-ddci*co*sO) - sf*(- ddci*so*sO));
	p.vy = v0*((e+cf)*(ddci*co*cO)  - sf*(ddci*so*cO));
	p.vz = v0*((e+cf)*co*ddsi - sf*ddsi*so);

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_dOmega(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double r = a*(1.-e*e)/(1. + e*cos(f));
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 

	double dcO = -sin(Omega);
	double dsO = cos(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	
	p.x = r*(dcO*(co*cf-so*sf) - dsO*(so*cf+co*sf)*ci);
	p.y = r*(dsO*(co*cf-so*sf) + dcO*(so*cf+co*sf)*ci);
	p.z = 0.;

	p.vx = v0*((e+cf)*(-ci*co*dsO - dcO*so) - sf*(co*dcO - ci*so*dsO));
	p.vy = v0*((e+cf)*(ci*co*dcO - dsO*so)  - sf*(co*dsO + ci*so*dcO));
	p.vz = 0.;

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_ddOmega(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double r = a*(1.-e*e)/(1. + e*cos(f));
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 

	double ddcO = -cos(Omega);
	double ddsO = -sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	
	p.x = r*(ddcO*(co*cf-so*sf) - ddsO*(so*cf+co*sf)*ci);
	p.y = r*(ddsO*(co*cf-so*sf) + ddcO*(so*cf+co*sf)*ci);
	p.z = 0.;

	p.vx = v0*((e+cf)*(-ci*co*ddsO - ddcO*so) - sf*(co*ddcO - ci*so*ddsO));
	p.vy = v0*((e+cf)*(ci*co*ddcO - ddsO*so)  - sf*(co*ddsO + ci*so*ddcO));
	p.vz = 0.;

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_domega(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double r = a*(1.-e*e)/(1. + e*cos(f));
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double dco = -sin(omega);
	double dso = cos(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.x = r*(cO*(dco*cf-dso*sf) - sO*(dso*cf+dco*sf)*ci);
	p.y = r*(sO*(dco*cf-dso*sf) + cO*(dso*cf+dco*sf)*ci);
	p.z = r*(dso*cf+dco*sf)*si;

	p.vx = v0*((e+cf)*(-ci*dco*sO - cO*dso) - sf*(dco*cO - ci*dso*sO));
	p.vy = v0*((e+cf)*(ci*dco*cO - sO*dso)  - sf*(dco*sO + ci*dso*cO));
	p.vz = v0*((e+cf)*dco*si - sf*si*dso);
	
	return p;
}

struct reb_particle reb_tools_orbit_to_particle_ddomega(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double r = a*(1.-e*e)/(1. + e*cos(f));
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double ddco = -cos(omega);
	double ddso = -sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.x = r*(cO*(ddco*cf-ddso*sf) - sO*(ddso*cf+ddco*sf)*ci);
	p.y = r*(sO*(ddco*cf-ddso*sf) + cO*(ddso*cf+ddco*sf)*ci);
	p.z = r*(ddso*cf+ddco*sf)*si;

	p.vx = v0*((e+cf)*(-ci*ddco*sO - cO*ddso) - sf*(ddco*cO - ci*ddso*sO));
	p.vy = v0*((e+cf)*(ci*ddco*cO - sO*ddso)  - sf*(ddco*sO + ci*ddso*cO));
	p.vz = v0*((e+cf)*ddco*si - sf*si*ddso);
	
	return p;
}

struct reb_particle reb_tools_orbit_to_particle_df(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double r = a*(1.-e*e)/(1. + e*cos(f));
	double dr = a*(1.-e*e)/((1. + e*cos(f))*(1. + e*cos(f)))*e*sin(f);
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double dcf = -sin(f);
	double dsf = cos(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.x = dr*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
	p.y = dr*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
	p.z = dr*(so*cf+co*sf)*si;
	
    p.x += r*(cO*(co*dcf-so*dsf) - sO*(so*dcf+co*dsf)*ci);
	p.y += r*(sO*(co*dcf-so*dsf) + cO*(so*dcf+co*dsf)*ci);
	p.z += r*(so*dcf+co*dsf)*si;

	p.vx = v0*(dcf*(-ci*co*sO - cO*so) - dsf*(co*cO - ci*so*sO));
	p.vy = v0*(dcf*(ci*co*cO - sO*so)  - dsf*(co*sO + ci*so*cO));
	p.vz = v0*(dcf*co*si - dsf*si*so);
	
	return p;
}

struct reb_particle reb_tools_orbit_to_particle_ddf(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double r = a*(1.-e*e)/(1. + e*cos(f));
	double dr = a*(1.-e*e)/((1. + e*cos(f))*(1. + e*cos(f)))*e*sin(f);
	double ddr = 2.*a*(1.-e*e)/((1. + e*cos(f))*(1. + e*cos(f))*(1. + e*cos(f)))*e*e*sin(f)*sin(f) + a*(1.-e*e)*e*cos(f)/((1. + e*cos(f))*(1. + e*cos(f)));
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double dcf = -sin(f);
	double dsf = cos(f);
	double ddcf = -cos(f);
	double ddsf = -sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.x = ddr*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
	p.y = ddr*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
	p.z = ddr*(so*cf+co*sf)*si;
    
    p.x += 2.*dr*(cO*(co*dcf-so*dsf) - sO*(so*dcf+co*dsf)*ci);
	p.y += 2.*dr*(sO*(co*dcf-so*dsf) + cO*(so*dcf+co*dsf)*ci);
	p.z += 2.*dr*(so*dcf+co*dsf)*si;
	
    p.x += r*(cO*(co*ddcf-so*ddsf) - sO*(so*ddcf+co*ddsf)*ci);
	p.y += r*(sO*(co*ddcf-so*ddsf) + cO*(so*ddcf+co*ddsf)*ci);
	p.z += r*(so*ddcf+co*ddsf)*si;

	p.vx = v0*(ddcf*(-ci*co*sO - cO*so) - ddsf*(co*cO - ci*so*sO));
	p.vy = v0*(ddcf*(ci*co*cO - sO*so)  - ddsf*(co*sO + ci*so*cO));
	p.vz = v0*(ddcf*co*si - ddsf*si*so);
	
	return p;
}

struct reb_particle reb_tools_orbit_to_particle_dm(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	p.m = 1.;
	double dv0 = 0.5*sqrt(G/a/(1.-e*e))/sqrt(m+primary.m); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);

	p.vx = dv0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
	p.vy = dv0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
	p.vz = dv0*((e+cf)*co*si - sf*si*so);

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_ddm(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double ddv0 = -0.25*sqrt(G/a/(1.-e*e))/sqrt(m+primary.m)/(m+primary.m); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);

	p.vx = ddv0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
	p.vy = ddv0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
	p.vz = ddv0*((e+cf)*co*si - sf*si*so);

	return p;
}


// Crossterms

struct reb_particle reb_tools_orbit_to_particle_da_dm(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){

	struct reb_particle p = {0};
	double ddv0 = -0.25/sqrt((m+primary.m)*a*a*a)*sqrt(G/(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.vx = ddv0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
	p.vy = ddv0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
	p.vz = ddv0*((e+cf)*co*si - sf*si*so);
	
	return p;
}



struct reb_particle reb_tools_orbit_to_particle_da_de(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
    double cosf = cos(f);
	double ddr = -(cosf*e*e+cosf+2.*e)/((cosf*e+1.)*(cosf*e+1.));
	double dv0_da = -0.5/sqrt(G*(m+primary.m)/a/(1.-e*e))*G*(m+primary.m)/(a*a)/(1.-e*e); 
	
	double dv0_da_de = e*dv0_da/(1.-e*e); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.x = ddr*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
	p.y = ddr*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
	p.z = ddr*(so*cf+co*sf)*si;

	p.vx = dv0_da_de*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
	p.vy = dv0_da_de*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
	p.vz = dv0_da_de*((e+cf)*co*si - sf*si*so);
	
    p.vx += dv0_da*(-ci*co*sO - cO*so);
	p.vy += dv0_da*(ci*co*cO - sO*so);
	p.vz += dv0_da*(co*si);

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_da_di(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){

	struct reb_particle p = {0};
	double dr = (1.-e*e)/(1. + e*cos(f));
	double dv0 = -0.5/sqrt(a*a*a)*sqrt(G*(m+primary.m)/(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double dci = -sin(inc);
	double dsi = cos(inc);
	
	p.x = dr*(- sO*(so*cf+co*sf)*dci);
	p.y = dr*(+ cO*(so*cf+co*sf)*dci);
	p.z = dr*(so*cf+co*sf)*dsi;

	p.vx = dv0*((e+cf)*(-dci*co*sO) - sf*(- dci*so*sO));
	p.vy = dv0*((e+cf)*(dci*co*cO)  - sf*(dci*so*cO));
	p.vz = dv0*((e+cf)*co*dsi - sf*dsi*so);
	
	return p;
}

struct reb_particle reb_tools_orbit_to_particle_da_dOmega(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){

	struct reb_particle p = {0};
	double dr = (1.-e*e)/(1. + e*cos(f));
	double dv0 = -0.5/sqrt(a*a*a)*sqrt(G*(m+primary.m)/(1.-e*e)); 

	double dcO = -sin(Omega);
	double dsO = cos(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	
	p.x = dr*(dcO*(co*cf-so*sf) - dsO*(so*cf+co*sf)*ci);
	p.y = dr*(dsO*(co*cf-so*sf) + dcO*(so*cf+co*sf)*ci);

	p.vx = dv0*((e+cf)*(-ci*co*dsO - dcO*so) - sf*(co*dcO - ci*so*dsO));
	p.vy = dv0*((e+cf)*(ci*co*dcO - dsO*so)  - sf*(co*dsO + ci*so*dcO));
	
	return p;
}

struct reb_particle reb_tools_orbit_to_particle_da_domega(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){

	struct reb_particle p = {0};
	double dr = (1.-e*e)/(1. + e*cos(f));
	double dv0 = -0.5/sqrt(a*a*a)*sqrt(G*(m+primary.m)/(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double dco = -sin(omega);
	double dso = cos(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.x = dr*(cO*(dco*cf-dso*sf) - sO*(dso*cf+dco*sf)*ci);
	p.y = dr*(sO*(dco*cf-dso*sf) + cO*(dso*cf+dco*sf)*ci);
	p.z = dr*(dso*cf+dco*sf)*si;

	p.vx = dv0*((e+cf)*(-ci*dco*sO - cO*dso) - sf*(dco*cO - ci*dso*sO));
	p.vy = dv0*((e+cf)*(ci*dco*cO - sO*dso)  - sf*(dco*sO + ci*dso*cO));
	p.vz = dv0*((e+cf)*dco*si - sf*si*dso);
	
	return p;
}
struct reb_particle reb_tools_orbit_to_particle_da_df(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){

	struct reb_particle p = {0};
	double dr = (1.-e*e)/(1. + e*cos(f));
	double ddr = e*sin(f)*(1.-e*e)/(1. + e*cos(f))/(1. + e*cos(f));
	double dv0 = -0.5/sqrt(a*a*a)*sqrt(G*(m+primary.m)/(1.-e*e));

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double dcf = -sin(f);
	double dsf = cos(f);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.x = dr*(cO*(co*dcf-so*dsf) - sO*(so*dcf+co*dsf)*ci);
	p.y = dr*(sO*(co*dcf-so*dsf) + cO*(so*dcf+co*dsf)*ci);
	p.z = dr*(so*dcf+co*dsf)*si;
	
    p.x += ddr*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
	p.y += ddr*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
	p.z += ddr*(so*cf+co*sf)*si;

	p.vx = dv0*(dcf*(-ci*co*sO - cO*so) - dsf*(co*cO - ci*so*sO));
	p.vy = dv0*(dcf*(ci*co*cO - sO*so)  - dsf*(co*sO + ci*so*cO));
	p.vz = dv0*(dcf*co*si - dsf*si*so);
	
	return p;
}
struct reb_particle reb_tools_orbit_to_particle_de_di(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
    double cosf = cos(f);
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 
	double dr = -a*(cosf*e*e+cosf+2.*e)/((cosf*e+1.)*(cosf*e+1.));
	double dv0 = sqrt(G*(m+primary.m)/a)*e/((1.-e*e)*sqrt(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double dci = -sin(inc);
	double dsi = cos(inc);
	
	p.x = dr*(- sO*(so*cf+co*sf)*dci);
	p.y = dr*(+ cO*(so*cf+co*sf)*dci);
	p.z = dr*(so*cf+co*sf)*dsi;

	p.vx = dv0*((e+cf)*(-dci*co*sO) - sf*(- dci*so*sO));
	p.vy = dv0*((e+cf)*(dci*co*cO)  - sf*(+ dci*so*cO));
	p.vz = dv0*((e+cf)*co*dsi - sf*dsi*so);
	
    p.vx += v0*(-dci*co*sO);
	p.vy += v0*(dci*co*cO);
	p.vz += v0*(co*dsi);

	return p;
}


struct reb_particle reb_tools_orbit_to_particle_de_dOmega(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
    double cosf = cos(f);
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 
	double dr = -a*(cosf*e*e+cosf+2.*e)/((cosf*e+1.)*(cosf*e+1.));
	double dv0 = sqrt(G*(m+primary.m)/a)*e/((1.-e*e)*sqrt(1.-e*e)); 

	double dcO = -sin(Omega);
	double dsO = cos(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	
	p.x = dr*(dcO*(co*cf-so*sf) - dsO*(so*cf+co*sf)*ci);
	p.y = dr*(dsO*(co*cf-so*sf) + dcO*(so*cf+co*sf)*ci);

	p.vx = dv0*((e+cf)*(-ci*co*dsO - dcO*so) - sf*(co*dcO - ci*so*dsO));
	p.vy = dv0*((e+cf)*(ci*co*dcO - dsO*so)  - sf*(co*dsO + ci*so*dcO));
	
    p.vx += v0*(-ci*co*dsO - dcO*so);
	p.vy += v0*(ci*co*dcO - dsO*so);

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_de_domega(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
    double cosf = cos(f);
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 
	double dr = -a*(cosf*e*e+cosf+2.*e)/((cosf*e+1.)*(cosf*e+1.));
	double dv0 = sqrt(G*(m+primary.m)/a)*e/((1.-e*e)*sqrt(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double dco = -sin(omega);
	double dso = cos(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.x = dr*(cO*(dco*cf-dso*sf) - sO*(dso*cf+dco*sf)*ci);
	p.y = dr*(sO*(dco*cf-dso*sf) + cO*(dso*cf+dco*sf)*ci);
	p.z = dr*(dso*cf+dco*sf)*si;

	p.vx = dv0*((e+cf)*(-ci*dco*sO - cO*dso) - sf*(dco*cO - ci*dso*sO));
	p.vy = dv0*((e+cf)*(ci*dco*cO - sO*dso)  - sf*(dco*sO + ci*dso*cO));
	p.vz = dv0*((e+cf)*dco*si - sf*si*dso);
	
    p.vx += v0*(-ci*dco*sO - cO*dso);
	p.vy += v0*(ci*dco*cO - sO*dso);
	p.vz += v0*(dco*si);

	return p;
}
struct reb_particle reb_tools_orbit_to_particle_de_df(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
    double cosf = cos(f);
	double dr = -a*(cosf*e*e+cosf+2.*e)/((cosf*e+1.)*(cosf*e+1.));
	double ddr = -a*(-sin(f)*e*e-sin(f))/((cosf*e+1.)*(cosf*e+1.))
	            -2.*e*sin(f) * a*(cosf*e*e+cosf+2.*e)/((cosf*e+1.)*(cosf*e+1.)*(cosf*e+1.));
	double dv0 = sqrt(G*(m+primary.m)/a)*e/((1.-e*e)*sqrt(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double dcf = -sin(f);
	double dsf = cos(f);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.x = dr*(cO*(co*dcf-so*dsf) - sO*(so*dcf+co*dsf)*ci);
	p.y = dr*(sO*(co*dcf-so*dsf) + cO*(so*dcf+co*dsf)*ci);
	p.z = dr*(so*dcf+co*dsf)*si;
	
    p.x += ddr*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
	p.y += ddr*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
	p.z += ddr*(so*cf+co*sf)*si;
	
    p.vx = dv0*(dcf*(-ci*co*sO - cO*so) - dsf*(co*cO - ci*so*sO));
	p.vy = dv0*(dcf*(ci*co*cO - sO*so)  - dsf*(co*sO + ci*so*cO));
	p.vz = dv0*(dcf*co*si - dsf*si*so);
	
	return p;
}
struct reb_particle reb_tools_orbit_to_particle_de_dm(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double dv0m = 0.5*G/a/(1.-e*e)/sqrt(G*(m+primary.m)/a/(1.-e*e)); 
	double dv0ea = 0.5*G/a/sqrt(G*(m+primary.m)/a)*e/((1.-e*e)*sqrt(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.vx = dv0ea*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
	p.vy = dv0ea*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
	p.vz = dv0ea*((e+cf)*co*si - sf*si*so);
	
    p.vx += dv0m*(-ci*co*sO - cO*so);
	p.vy += dv0m*(ci*co*cO - sO*so);
	p.vz += dv0m*(co*si);

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_di_dOmega(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double r = a*(1.-e*e)/(1. + e*cos(f));
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 

	double dcO = -sin(Omega);
	double dsO = cos(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double dci = -sin(inc);
	
	p.x = r*(- dsO*(so*cf+co*sf)*dci);
	p.y = r*(+ dcO*(so*cf+co*sf)*dci);

	p.vx = v0*((e+cf)*(-dci*co*dsO) - sf*(- dci*so*dsO));
	p.vy = v0*((e+cf)*(dci*co*dcO)  - sf*(dci*so*dcO));

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_di_domega(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double r = a*(1.-e*e)/(1. + e*cos(f));
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double dco = -sin(omega);
	double dso = cos(omega);
	double cf = cos(f);
	double sf = sin(f);
	double dci = -sin(inc);
	double dsi = cos(inc);
	
	p.x = r*(- sO*(dso*cf+dco*sf)*dci);
	p.y = r*(+ cO*(dso*cf+dco*sf)*dci);
	p.z = r*(dso*cf+dco*sf)*dsi;

	p.vx = v0*((e+cf)*(-dci*dco*sO) - sf*(- dci*dso*sO));
	p.vy = v0*((e+cf)*(dci*dco*cO)  - sf*(dci*dso*cO));
	p.vz = v0*((e+cf)*dco*dsi - sf*dsi*dso);

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_di_df(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double r = a*(1.-e*e)/(1. + e*cos(f));
	double dr = e*sin(f)*a*(1.-e*e)/(1. + e*cos(f))/(1. + e*cos(f));
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double dcf = -sin(f);
	double dsf = cos(f);
	double cf = cos(f);
	double sf = sin(f);
	double dci = -sin(inc);
	double dsi = cos(inc);
	
	p.x = r*(- sO*(so*dcf+co*dsf)*dci);
	p.y = r*(+ cO*(so*dcf+co*dsf)*dci);
	p.z = r*(so*dcf+co*dsf)*dsi;
	
    p.x += dr*(- sO*(so*cf+co*sf)*dci);
	p.y += dr*(+ cO*(so*cf+co*sf)*dci);
	p.z += dr*(so*cf+co*sf)*dsi;

	p.vx = v0*(dcf*(-dci*co*sO) - dsf*(- dci*so*sO));
	p.vy = v0*(dcf*(dci*co*cO)  - dsf*(dci*so*cO));
	p.vz = v0*(dcf*co*dsi - dsf*dsi*so);

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_di_dm(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double dv0 = 0.5/sqrt(m+primary.m)*sqrt(G/a/(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double dci = -sin(inc);
	double dsi = cos(inc);
	
	p.vx = dv0*((e+cf)*(-dci*co*sO) - sf*(- dci*so*sO));
	p.vy = dv0*((e+cf)*(dci*co*cO)  - sf*(dci*so*cO));
	p.vz = dv0*((e+cf)*co*dsi - sf*dsi*so);
	
	return p;
}

struct reb_particle reb_tools_orbit_to_particle_dOmega_domega(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double r = a*(1.-e*e)/(1. + e*cos(f));
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 

	double dcO = -sin(Omega);
	double dsO = cos(Omega);
	double dco = -sin(omega);
	double dso = cos(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	
	p.x = r*(dcO*(dco*cf-dso*sf) - dsO*(dso*cf+dco*sf)*ci);
	p.y = r*(dsO*(dco*cf-dso*sf) + dcO*(dso*cf+dco*sf)*ci);

	p.vx = v0*((e+cf)*(-ci*dco*dsO - dcO*dso) - sf*(dco*dcO - ci*dso*dsO));
	p.vy = v0*((e+cf)*(ci*dco*dcO - dsO*dso)  - sf*(dco*dsO + ci*dso*dcO));

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_dOmega_df(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double r = a*(1.-e*e)/(1. + e*cos(f));
	double dr = e*sin(f)*a*(1.-e*e)/(1. + e*cos(f))/(1. + e*cos(f));
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 

	double dcO = -sin(Omega);
	double dsO = cos(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double dcf = -sin(f);
	double dsf = cos(f);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	
	p.x = r*(dcO*(co*dcf-so*dsf) - dsO*(so*dcf+co*dsf)*ci);
	p.y = r*(dsO*(co*dcf-so*dsf) + dcO*(so*dcf+co*dsf)*ci);
	
    p.x += dr*(dcO*(co*cf-so*sf) - dsO*(so*cf+co*sf)*ci);
	p.y += dr*(dsO*(co*cf-so*sf) + dcO*(so*cf+co*sf)*ci);

	p.vx = v0*((dcf)*(-ci*co*dsO - dcO*so) - dsf*(co*dcO - ci*so*dsO));
	p.vy = v0*((dcf)*(ci*co*dcO - dsO*so)  - dsf*(co*dsO + ci*so*dcO));

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_dOmega_dm(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double dv0 = 0.5/sqrt(m+primary.m)*sqrt(G/a/(1.-e*e)); 

	double dcO = -sin(Omega);
	double dsO = cos(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	
	p.vx = dv0*((e+cf)*(-ci*co*dsO - dcO*so) - sf*(co*dcO - ci*so*dsO));
	p.vy = dv0*((e+cf)*(ci*co*dcO - dsO*so)  - sf*(co*dsO + ci*so*dcO));

	return p;
}

struct reb_particle reb_tools_orbit_to_particle_domega_df(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double r = a*(1.-e*e)/(1. + e*cos(f));
	double dr = e*sin(f)*a*(1.-e*e)/(1. + e*cos(f))/(1. + e*cos(f));
	double v0 = sqrt(G*(m+primary.m)/a/(1.-e*e)); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double dco = -sin(omega);
	double dso = cos(omega);
	double dcf = -sin(f);
	double dsf = cos(f);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.x = r*(cO*(dco*dcf-dso*dsf) - sO*(dso*dcf+dco*dsf)*ci);
	p.y = r*(sO*(dco*dcf-dso*dsf) + cO*(dso*dcf+dco*dsf)*ci);
	p.z = r*(dso*dcf+dco*dsf)*si;
	
    p.x += dr*(cO*(dco*cf-dso*sf) - sO*(dso*cf+dco*sf)*ci);
	p.y += dr*(sO*(dco*cf-dso*sf) + cO*(dso*cf+dco*sf)*ci);
	p.z += dr*(dso*cf+dco*sf)*si;

	p.vx = v0*((dcf)*(-ci*dco*sO - cO*dso) - dsf*(dco*cO - ci*dso*sO));
	p.vy = v0*((dcf)*(ci*dco*cO - sO*dso)  - dsf*(dco*sO + ci*dso*cO));
	p.vz = v0*((dcf)*dco*si - dsf*si*dso);
	
	return p;
}

struct reb_particle reb_tools_orbit_to_particle_domega_dm(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double dv0 = 0.5*sqrt(G/a/(1.-e*e))/sqrt(m+primary.m); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double dco = -sin(omega);
	double dso = cos(omega);
	double cf = cos(f);
	double sf = sin(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.vx = dv0*((e+cf)*(-ci*dco*sO - cO*dso) - sf*(dco*cO - ci*dso*sO));
	p.vy = dv0*((e+cf)*(ci*dco*cO - sO*dso)  - sf*(dco*sO + ci*dso*cO));
	p.vz = dv0*((e+cf)*dco*si - sf*si*dso);
	
	return p;
}

struct reb_particle reb_tools_orbit_to_particle_df_dm(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	struct reb_particle p = {0};
	double dv0 = 0.5*sqrt(G/a/(1.-e*e))/sqrt(m+primary.m); 

	double cO = cos(Omega);
	double sO = sin(Omega);
	double co = cos(omega);
	double so = sin(omega);
	double dcf = -sin(f);
	double dsf = cos(f);
	double ci = cos(inc);
	double si = sin(inc);
	
	p.vx = dv0*(dcf*(-ci*co*sO - cO*so) - dsf*(co*cO - ci*so*sO));
	p.vy = dv0*(dcf*(ci*co*cO - sO*so)  - dsf*(co*sO + ci*so*cO));
	p.vz = dv0*(dcf*co*si - dsf*si*so);
	
	return p;
}
