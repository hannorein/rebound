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
#include <string.h>
#include <sys/time.h>
#include <stdint.h>
#include "particle.h"
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
    const int N_var = r->N_var;
    const int _N_active = ((r->N_active==-1)?N:r->N_active) - N_var;
    const struct reb_particle* restrict const particles = r->particles;
    double e_kin = 0.;
    double e_pot = 0.;
    int N_interact = (r->testparticle_type==0)?_N_active:(N-N_var);
    for (int i=0;i<N_interact;i++){
        struct reb_particle pi = particles[i];
        e_kin += 0.5 * pi.m * (pi.vx*pi.vx + pi.vy*pi.vy + pi.vz*pi.vz);
    }
    for (int i=0;i<_N_active;i++){
        struct reb_particle pi = particles[i];
        for (int j=i+1;j<N_interact;j++){
            struct reb_particle pj = particles[j];
            double dx = pi.x - pj.x;
            double dy = pi.y - pj.y;
            double dz = pi.z - pj.z;
            e_pot -= r->G*pj.m*pi.m/sqrt(dx*dx + dy*dy + dz*dz);
        }
    }
    
    return e_kin + e_pot + r->energy_offset;
}

struct reb_vec3d reb_tools_angular_momentum(const struct reb_simulation* const r){
	const int N = r->N;
	const struct reb_particle* restrict const particles = r->particles;
	const int N_var = r->N_var;
    struct reb_vec3d L = {0};
    for (int i=0;i<N-N_var;i++){
		struct reb_particle pi = particles[i];
        L.x += pi.m*(pi.y*pi.vz - pi.z*pi.vy);
        L.y += pi.m*(pi.z*pi.vx - pi.x*pi.vz);
        L.z += pi.m*(pi.x*pi.vy - pi.y*pi.vx);
	}
	return L;
}

void reb_move_to_hel(struct reb_simulation* const r){
    const int N_real = r->N - r->N_var;
    if (N_real>0){
	    struct reb_particle* restrict const particles = r->particles;
        struct reb_particle hel = r->particles[0];
        // Note: Variational particles will not be affected.
        for (int i=1;i<N_real;i++){
            particles[i].x  -= hel.x;
            particles[i].y  -= hel.y;
            particles[i].z  -= hel.z;
            particles[i].vx -= hel.vx;
            particles[i].vy -= hel.vy;
            particles[i].vz -= hel.vz;
        }
        r->particles[0].x = 0.;
        r->particles[0].y = 0.;
        r->particles[0].z = 0.;
        r->particles[0].vx = 0.;
        r->particles[0].vy = 0.;
        r->particles[0].vz = 0.;
    }
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

void reb_serialize_particle_data(struct reb_simulation* r, uint32_t* hash, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]){
    const int N_real = r->N - r->N_var;
    struct reb_particle* restrict const particles = r->particles;
    for (int i=0;i<N_real;i++){
        if (hash){
            hash[i] = particles[i].hash;
        }
        if (m){
            m[i] = particles[i].m;
        }
        if (radius){
            radius[i] = particles[i].r;
        }
        if (xyz){
            xyz[i][0] = particles[i].x;
            xyz[i][1] = particles[i].y;
            xyz[i][2] = particles[i].z;
        }
        if (vxvyvz){
            vxvyvz[i][0] = particles[i].vx;
            vxvyvz[i][1] = particles[i].vy;
            vxvyvz[i][2] = particles[i].vz;
        }
        if (xyzvxvyvz){
            xyzvxvyvz[i][0] = particles[i].x;
            xyzvxvyvz[i][1] = particles[i].y;
            xyzvxvyvz[i][2] = particles[i].z;
            xyzvxvyvz[i][3] = particles[i].vx;
            xyzvxvyvz[i][4] = particles[i].vy;
            xyzvxvyvz[i][5] = particles[i].vz;
        }
    }
}

void reb_set_serialized_particle_data(struct reb_simulation* r, uint32_t* hash, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]){
    const int N_real = r->N - r->N_var;
    struct reb_particle* restrict const particles = r->particles;
    for (int i=0;i<N_real;i++){
        if (hash){
           particles[i].hash = hash[i];
        }
        if (m){
            particles[i].m = m[i];
        }
        if (radius){
            particles[i].r = radius[i] ;
        }
        if (xyz){
            particles[i].x = xyz[i][0];
            particles[i].y = xyz[i][1];
            particles[i].z = xyz[i][2];
        }
        if (vxvyvz){
            particles[i].vx = vxvyvz[i][0];
            particles[i].vy = vxvyvz[i][1];
            particles[i].vz = vxvyvz[i][2];
        }
        if (xyzvxvyvz){
            particles[i].x = xyzvxvyvz[i][0];
            particles[i].y = xyzvxvyvz[i][1];
            particles[i].z = xyzvxvyvz[i][2];
            particles[i].vx = xyzvxvyvz[i][3];
            particles[i].vy = xyzvxvyvz[i][4];
            particles[i].vz = xyzvxvyvz[i][5];
        }
    }
}

struct reb_particle reb_get_com_of_pair(struct reb_particle p1, struct reb_particle p2){
	p1.x   = p1.x*p1.m + p2.x*p2.m;		
	p1.y   = p1.y*p1.m + p2.y*p2.m;
	p1.z   = p1.z*p1.m + p2.z*p2.m;
	p1.vx  = p1.vx*p1.m + p2.vx*p2.m;
	p1.vy  = p1.vy*p1.m + p2.vy*p2.m;
	p1.vz  = p1.vz*p1.m + p2.vz*p2.m;
	p1.ax  = p1.ax*p1.m + p2.ax*p2.m;
	p1.ay  = p1.ay*p1.m + p2.ay*p2.m;
	p1.az  = p1.az*p1.m + p2.az*p2.m;
    
	p1.m  += p2.m;
	if (p1.m>0.){
		p1.x  /= p1.m;
		p1.y  /= p1.m;
		p1.z  /= p1.m;
		p1.vx /= p1.m;
		p1.vy /= p1.m;
		p1.vz /= p1.m;
		p1.ax /= p1.m;
		p1.ay /= p1.m;
		p1.az /= p1.m;
	}
	return p1;
}

struct reb_particle reb_get_com_without_particle(struct reb_particle com, struct reb_particle p){
    com.x = com.x*com.m - p.x*p.m;
    com.y = com.y*com.m - p.y*p.m;
    com.z = com.z*com.m - p.z*p.m;
    com.vx = com.vx*com.m - p.vx*p.m;
    com.vy = com.vy*com.m - p.vy*p.m;
    com.vz = com.vz*com.m - p.vz*p.m;
    com.ax = com.ax*com.m - p.ax*p.m;
    com.ay = com.ay*com.m - p.ay*p.m;
    com.az = com.az*com.m - p.az*p.m;
    com.m -= p.m; 

    if (com.m > 0.){
        com.x /= com.m;
        com.y /= com.m;
        com.z /= com.m;
        com.vx /= com.m;
        com.vy /= com.m;
        com.vz /= com.m;
        com.ax /= com.m;
        com.ay /= com.m;
        com.az /= com.m;
    }
    return com;
}

struct reb_particle reb_get_com_range(struct reb_simulation* r, int first, int last){
	struct reb_particle com = {0};
	for(int i=first; i<last; i++){
		com = reb_get_com_of_pair(com, r->particles[i]);
	}
	return com;
}

struct reb_particle reb_get_com(struct reb_simulation* r){
    int N_real = r->N-r->N_var;
	return reb_get_com_range(r, 0, N_real); 
}

struct reb_particle reb_get_jacobi_com(struct reb_particle* p){
	int p_index = reb_get_particle_index(p);
	struct reb_simulation* r = p->sim;
    return reb_get_com_range(r, 0, p_index);
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

static double mod2pi(double f){
	while(f < 0.){
		f += 2.*M_PI;
	}
	while(f > 2.*M_PI){
		f -= 2.*M_PI;
	}
	return f;
}

double reb_tools_M_to_E(double e, double M){
	double E;
	if(e < 1.){
        M = mod2pi(M); // avoid numerical artefacts for negative numbers
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
	double E = reb_tools_M_to_E(e, M);
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

#define TINY 1.E-308 		///< Close to smallest representable floating point number, used for orbit calculation

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
    if(primary.m < TINY){
        *err = 6;       // Primary has no mass.
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
    o.rhill = nan("");

    return o;
}

#define MIN_REL_ERROR 1.0e-12	///< Close to smallest relative floating point number, used for orbit calculation
#define MIN_INC 1.e-8		///< Below this inclination, the broken angles pomega and theta equal the corresponding 
							///< unbroken angles to within machine precision, so a practical boundary for planar orbits
							//
#define MIN_ECC 1.e-8       ///< Below this eccentricity, corrections at order e^2 are below machine precision, so we use
                            ///< stable expressions accurate to O(e) for the mean longitude below for near-circular orbits.
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
    
    o.rhill = o.a*cbrt(p.m/(3.*primary.m));
	
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
		o.theta = acos2(dx, o.d, dy);		// cos theta is dot product of x and r vectors (true longitude). 
        o.pomega = acos2(ex, o.e, ey);		// cos pomega is dot product of x and e unit vectors.  Will = 0 if e=0.

		if(o.inc < M_PI/2.){
			o.omega = o.pomega - o.Omega;
			o.f = o.theta - o.pomega;
            if(o.e > MIN_ECC){              // pomega well defined
			    o.l = o.pomega + o.M;
            }
            else{                           // when e << 1 and pomega ill defined, use l = theta+(M-f). M-f is O(e) so well behaved
                o.l = o.theta - 2.*o.e*sin(o.f); // M-f from Murray & Dermott Eq 2.93. This way l->theta smoothly as e->0
            }
		}
		else{
			o.omega = o.Omega - o.pomega;
			o.f = o.pomega - o.theta;
            if(o.e > MIN_ECC){              // pomega well defined
			    o.l = o.pomega - o.M;
            }
            else{                           // when e << 1 and pomega ill defined, use l = theta+(M-f). M-f is O(e) so well behaved
                o.l = o.theta + 2.*o.e*sin(o.f); // M-f from Murray & Dermott Eq 2.93 (retrograde changes sign). This way l->theta smoothly as e->0
            }
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
            if(o.e > MIN_ECC){              // pomega well defined
			    o.l = o.pomega + o.M;
            }
            else{                           // when e << 1 and pomega ill defined, use l = theta+(M-f). M-f is O(e) so well behaved
                o.l = o.theta - 2.*o.e*sin(o.f); // M-f from Murray & Dermott Eq 2.93. This way l->theta smoothly as e->0
            }
		}
		else{
			o.pomega = o.Omega - o.omega;
			o.f = wpf - o.omega;
			o.theta = o.Omega - wpf;
            if(o.e > MIN_ECC){              // pomega well defined
			    o.l = o.pomega - o.M;
            }
            else{                           // when e << 1 and pomega ill defined, use l = theta+(M-f). M-f is O(e) so well behaved
                o.l = o.theta + 2.*o.e*sin(o.f); // M-f from Murray & Dermott Eq 2.93 (retrograde changes sign). This way l->theta smoothly as e->0
            }
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


void reb_tools_solve_kepler_pal(double h, double k, double lambda, double* p, double* q){
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
        double E = reb_tools_M_to_E(e, M);
        *p = e*sin(E); 
        *q = e*cos(E); 
    }
}

void reb_tools_particle_to_pal(double G, struct reb_particle p, struct reb_particle primary, double *a, double* lambda, double* k, double* h, double* ix, double* iy){
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

struct reb_particle reb_tools_pal_to_particle(double G, struct reb_particle primary, double m, double a, double lambda, double k, double h, double ix, double iy){
    struct reb_particle np = {0.};
    np.m = m;

    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double xi = a*(clp + p/(2.-l)*h -k);
    double eta = a*(slp - p/(2.-l)*k -h);

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

#define ROT32(x, y) ((x << y) | (x >> (32 - y))) // avoid effort
static uint32_t reb_murmur3_32(const char *key, uint32_t len, uint32_t seed) {
    // Source: Wikipedia
    static const uint32_t c1 = 0xcc9e2d51;
    static const uint32_t c2 = 0x1b873593;
    static const uint32_t r1 = 15;
    static const uint32_t r2 = 13;
    static const uint32_t m = 5;
    static const uint32_t n = 0xe6546b64;

    uint32_t hash = seed;

    const int nblocks = len / 4;
    const uint32_t *blocks = (const uint32_t *) key;
    int i;
    uint32_t k;
    for (i = 0; i < nblocks; i++) {
        k = blocks[i];
        k *= c1;
        k = ROT32(k, r1);
        k *= c2;

        hash ^= k;
        hash = ROT32(hash, r2) * m + n;
    }

    const uint8_t *tail = (const uint8_t *) (key + nblocks * 4);
    uint32_t k1 = 0;

    switch (len & 3) {
    case 3:
        k1 ^= tail[2] << 16;
    case 2:
        k1 ^= tail[1] << 8;
    case 1:
        k1 ^= tail[0];

        k1 *= c1;
        k1 = ROT32(k1, r1);
        k1 *= c2;
        hash ^= k1;
    }

    hash ^= len;
    hash ^= (hash >> 16);
    hash *= 0x85ebca6b;
    hash ^= (hash >> 13);
    hash *= 0xc2b2ae35;
    hash ^= (hash >> 16);

    return hash;
}

uint32_t reb_hash(const char* str){
    const int reb_seed = 1983;
    return reb_murmur3_32(str,(uint32_t)strlen(str),reb_seed);
}

void reb_simulation_imul(struct reb_simulation* r, double scalar_pos, double scalar_vel){
    const int N = r->N;
    struct reb_particle* restrict const particles = r->particles;
	for (int i=0;i<N;i++){
        particles[i].x *= scalar_pos;
        particles[i].y *= scalar_pos;
        particles[i].z *= scalar_pos;
        particles[i].vx *= scalar_vel;
        particles[i].vy *= scalar_vel;
        particles[i].vz *= scalar_vel;
    }
}

int reb_simulation_iadd(struct reb_simulation* r, struct reb_simulation* r2){
    const int N = r->N;
    const int N2 = r2->N;
    if (N!=N2) return -1;
    struct reb_particle* restrict const particles = r->particles;
    const struct reb_particle* restrict const particles2 = r2->particles;
	for (int i=0;i<N;i++){
        particles[i].x += particles2[i].x;
        particles[i].y += particles2[i].y;
        particles[i].z += particles2[i].z;
        particles[i].vx += particles2[i].vx;
        particles[i].vy += particles2[i].vy;
        particles[i].vz += particles2[i].vz;
    }
    return 0;
}

int reb_simulation_isub(struct reb_simulation* r, struct reb_simulation* r2){
    const int N = r->N;
    const int N2 = r2->N;
    if (N!=N2) return -1;
    struct reb_particle* restrict const particles = r->particles;
    const struct reb_particle* restrict const particles2 = r2->particles;
	for (int i=0;i<N;i++){
        particles[i].x -= particles2[i].x;
        particles[i].y -= particles2[i].y;
        particles[i].z -= particles2[i].z;
        particles[i].vx -= particles2[i].vx;
        particles[i].vy -= particles2[i].vy;
        particles[i].vz -= particles2[i].vz;
    }
    return 0;
}

