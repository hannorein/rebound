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
#include <string.h>
#ifdef _WIN32
#define strtok_r strtok_s
#define REB_RAND_MAX 2147483647  // INT_MAX
#else // Linux and MacOS
#define REB_RAND_MAX RAND_MAX
#endif // _WIN32
#include <stdint.h>
#include <stdarg.h>
#include "rebound.h"
#include "particle.h"
#include "rebound.h"
#include "tools.h"
#include "tree.h"
#include "boundary.h"
#ifdef MPI
#include "communication_mpi.h"
#endif // MPI
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b


unsigned int reb_tools_get_rand_seed(){
	struct reb_timeval tim;
	gettimeofday(&tim, NULL);
	return tim.tv_usec + getpid();
}

double reb_random_uniform(struct reb_simulation* r, double min, double max){
    unsigned int seed;
    unsigned int* seedp;
    if (r){
        seedp = &r->rand_seed;
    }else{
        seed = reb_tools_get_rand_seed();
        seedp = &seed;
    }
	return ((double)rand_r(seedp))/((double)(REB_RAND_MAX))*(max-min)+min;
}


double reb_random_powerlaw(struct reb_simulation* r, double min, double max, double slope){
	double y = reb_random_uniform(r, 0., 1.);
	if(slope == -1) return exp(y*log(max/min) + log(min));
    else return pow( (pow(max,slope+1.)-pow(min,slope+1.))*y+pow(min,slope+1.), 1./(slope+1.));
}

double reb_random_normal(struct reb_simulation* r, double variance){
	double v1=0.,v2=0.,rsq=1.;
    unsigned int seed;
    unsigned int* seedp;
    if (r){
        seedp = &r->rand_seed;
    }else{
        seed = reb_tools_get_rand_seed();
        seedp = &seed;
    }
	while(rsq>=1. || rsq<1.0e-12){
		v1=2.*((double)rand_r(seedp))/((double)(REB_RAND_MAX))-1.0;
		v2=2.*((double)rand_r(seedp))/((double)(REB_RAND_MAX))-1.0;
		rsq=v1*v1+v2*v2;
	}
	// Note: This gives another random variable for free, but we'll throw it away for simplicity and for thread-safety.
	return 	v1*sqrt(-2.*log(rsq)/rsq*variance);
}

double reb_random_rayleigh(struct reb_simulation* r, double sigma){
	double y = reb_random_uniform(r, 0.,1.);
	return sigma*sqrt(-2*log(y));
}

/// Other helper routines
double reb_simulation_energy(const struct reb_simulation* const r){
    const int N = r->N;
    const int N_var = r->N_var;
    const int _N_active = (r->N_active==-1)?(N-N_var):r->N_active;
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

struct reb_vec3d reb_simulation_angular_momentum(const struct reb_simulation* const r){
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

void reb_simulation_move_to_hel(struct reb_simulation* const r){
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


void reb_simulation_move_to_com(struct reb_simulation* const r){
	struct reb_particle com = reb_simulation_com(r); // Particles will be redistributed in this call if MPI used
	struct reb_particle* restrict const particles = r->particles;
    const int N_real = r->N - r->N_var; 
    
    // First do second order
    for (int v=0;v<r->N_var_config;v++){
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
    for (int v=0;v<r->N_var_config;v++){
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
    
    // Check boundaries and update tree if needed
    reb_boundary_check(r);     
    if (r->gravity==REB_GRAVITY_TREE || r->collision==REB_COLLISION_TREE || r->collision==REB_COLLISION_LINETREE){
        reb_simulation_update_tree(r);          
    }
#ifdef MPI
    reb_communication_mpi_distribute_particles(r);
#endif // MPI
}

void reb_simulation_get_serialized_particle_data(struct reb_simulation* r, uint32_t* hash, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]){
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

void reb_simulation_set_serialized_particle_data(struct reb_simulation* r, uint32_t* hash, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]){
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

struct reb_particle reb_particle_com_of_pair(struct reb_particle p1, struct reb_particle p2){
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

struct reb_particle reb_simulation_com_range(struct reb_simulation* r, int first, int last){
	struct reb_particle com = {0};
	for(int i=first; i<last; i++){
		com = reb_particle_com_of_pair(com, r->particles[i]);
	}
	return com;
}

struct reb_particle reb_simulation_com(struct reb_simulation* r){
#ifdef MPI
    reb_communication_mpi_distribute_particles(r);
    int N_real = r->N-r->N_var;
    struct reb_particle com_local = reb_simulation_com_range(r, 0, N_real);
	struct reb_particle com = {0};
    com_local.x  *= com_local.m;
    com_local.y  *= com_local.m;
    com_local.z  *= com_local.m;
    com_local.vx *= com_local.m;
    com_local.vy *= com_local.m;
    com_local.vz *= com_local.m;
    com_local.ax *= com_local.m;
    com_local.ay *= com_local.m;
    com_local.az *= com_local.m;

    MPI_Allreduce(&com_local.x, &com.x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&com_local.x, &com.y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&com_local.x, &com.z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&com_local.vx, &com.vx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&com_local.vx, &com.vy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&com_local.vx, &com.vz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&com_local.ax, &com.ax, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&com_local.ax, &com.ay, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&com_local.ax, &com.az, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&com_local.m, &com.m, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    if (com.m > 0){
        com.x  /= com.m;
        com.y  /= com.m;
        com.z  /= com.m;
        com.vx /= com.m;
        com.vy /= com.m;
        com.vz /= com.m;
        com.ax /= com.m;
        com.ay /= com.m;
        com.az /= com.m;
    }

	return com; 
#else // MPI
    int N_real = r->N-r->N_var;
    return reb_simulation_com_range(r, 0, N_real);
#endif // MPI
}

struct reb_particle reb_simulation_jacobi_com(struct reb_particle* p){
	int p_index = reb_simulation_particle_index(p);
	struct reb_simulation* r = p->sim;
    return reb_simulation_com_range(r, 0, p_index);
}
	
void reb_simulation_add_plummer(struct reb_simulation* r, int _N, double M, double R) {
	// Algorithm from:	
	// http://adsabs.harvard.edu/abs/1974A%26A....37..183A
	
	double E = 3./64.*M_PI*M*M/R;
	for (int i=0;i<_N;i++){
		struct reb_particle star = {0};
		double _r = pow(pow(reb_random_uniform(r, 0,1),-2./3.)-1.,-1./2.);
		double x2 = reb_random_uniform(r, 0,1);
		double x3 = reb_random_uniform(r, 0,2.*M_PI);
		star.z = (1.-2.*x2)*_r;
		star.x = sqrt(_r*_r-star.z*star.z)*cos(x3);
		star.y = sqrt(_r*_r-star.z*star.z)*sin(x3);
		double x5,g,q;
		do{
			x5 = reb_random_uniform(r, 0.,1.);
			q = reb_random_uniform(r, 0.,1.);
			g = q*q*pow(1.-q*q,7./2.);
		}while(0.1*x5>g);
		double ve = pow(2.,1./2.)*pow(1.+_r*_r,-1./4.);
		double v = q*ve;
		double x6 = reb_random_uniform(r, 0.,1.);
		double x7 = reb_random_uniform(r, 0.,2.*M_PI);
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

		reb_simulation_add(r, star);
	}
}

double reb_mod2pi(double f){
    const double pi2 = 2.*M_PI;
    return fmod(pi2 + fmod(f, pi2), pi2);
}

double reb_M_to_E(double e, double M){
	double E;
	if(e < 1.){
        M = reb_mod2pi(M); // avoid numerical artefacts for negative numbers
		E = e < 0.8 ? M : M_PI;
		double F = E - e*sin(E) - M;
		for(int i=0; i<100; i++){
			E = E - F/(1.-e*cos(E));
			F = E - e*sin(E) - M;
			if(fabs(F) < 1.e-16){
				break;
			}
		}
		E = reb_mod2pi(E);
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

double reb_E_to_f(double e, double E){
	if(e > 1.){
		return reb_mod2pi(2.*atan(sqrt((1.+e)/(e-1.))*tanh(0.5*E)));
	}
	else{
		return reb_mod2pi(2.*atan(sqrt((1.+e)/(1.-e))*tan(0.5*E)));
	}
}

double reb_M_to_f(double e, double M){
	double E = reb_M_to_E(e, M);
    return reb_E_to_f(e, E);
}

static const char* reb_string_for_particle_error(int err){
    if (err==1)
        return "Cannot set e exactly to 1.";
    if (err==2)
        return "Eccentricity must be greater than or equal to zero.";
    if (err==3)
        return "Bound orbit (a > 0) must have e < 1.";
    if (err==4)
        return "Unbound orbit (a < 0) must have e > 1.";
    if (err==5)
        return "Unbound orbit can't have f beyond the range allowed by the asymptotes set by the hyperbola.";
    if (err==6)
        return "Primary has no mass.";
    if (err==7)
        return "Cannot mix Pal coordinates (h,k,ix,iy) with certain orbital elements (e, inc, Omega, omega, pomega, f, M, E, theta, T). Use longitude l to indicate the phase.";
    if (err==8)
        return "Cannot pass cartesian coordinates and orbital elements (incl primary) at the same time.";
    if (err==9)
        return "Need to pass reb_simulation object when initializing particle with orbital elements.";
    if (err==10)
        return "Need to pass either semi-major axis or orbital period to initialize particle using orbital elements.";
    if (err==11)
        return "Need to pass either semi-major axis or orbital period, but not both.";
    if (err==12)
        return "(ix, iy) coordinates are not valid. Squared sum exceeds 4.";
    if (err==13)
        return "Cannot pass both (omega, pomega) together.";
    if (err==14)
        return "Can only pass one longitude/anomaly in the set (f, M, E, l, theta, T).";
    return "An unknown error occured during reb_simulation_add_fmt().";

}

static struct reb_particle reb_particle_from_fmt_errV(struct reb_simulation* r, int* err, const char* fmt, va_list args);

void reb_simulation_add_fmt(struct reb_simulation* r, const char* fmt, ...){
    if (!r){
        fprintf(stderr, "\n\033[1mError!\033[0m Simulation can't be NULL1.\n");
        return;
    }

    int err = 0;
    va_list args;
    va_start(args, fmt);
    struct reb_particle particle = reb_particle_from_fmt_errV(r, &err, fmt, args);
    va_end(args);

    if (err==0){ // Success
        reb_simulation_add(r, particle);
    }else{
        const char* error_string = reb_string_for_particle_error(err);
        reb_simulation_error(r, error_string);
    }
}

struct reb_particle reb_particle_from_fmt(struct reb_simulation* r, const char* fmt, ...){
    int err = 0;

    va_list args;
    va_start(args, fmt);
    struct reb_particle particle = reb_particle_from_fmt_errV(r, &err, fmt, args);
    va_end(args);
    
    if (err==0){ // Success
        return particle;
    }else{
        const char* error_string = reb_string_for_particle_error(err);
        fprintf(stderr, "\n\033[1mError!\033[0m %s\n", error_string);
        return reb_particle_nan();
    }
}

static struct reb_particle reb_particle_from_fmt_errV(struct reb_simulation* r, int* err, const char* fmt, va_list args){
    double m = 0;
    double radius = 0;
    uint32_t hash = 0;
    double x = nan("");
    double y = nan("");
    double z = nan("");
    double vx = nan("");
    double vy = nan("");
    double vz = nan("");
    double a = nan("");
    double P = nan("");
    double e = nan("");
    double inc = nan("");
    double Omega = nan("");
    double omega = nan("");
    double pomega = nan("");
    double f = nan("");
    double M = nan("");
    double E = nan("");
    double l = nan("");
    double theta = nan("");
    double T = nan("");
    double h = nan("");
    double k = nan("");
    double ix = nan("");
    double iy = nan("");
    struct reb_particle primary = {0};
    int primary_given = 0;

    char *sep = " \t\n,;";

    char* fmt_c = strdup(fmt);
    char* token;
    char* rest = fmt_c;

    while ((token = strtok_r(rest, sep, &rest))){
        if (0==strcmp(token,"m"))
            m = va_arg(args, double);
        if (0==strcmp(token,"r"))
            radius = va_arg(args, double);
        if (0==strcmp(token,"x"))
            x = va_arg(args, double);
        if (0==strcmp(token,"y"))
            y = va_arg(args, double);
        if (0==strcmp(token,"z"))
            z = va_arg(args, double);
        if (0==strcmp(token,"vx"))
            vx = va_arg(args, double);
        if (0==strcmp(token,"vy"))
            vy = va_arg(args, double);
        if (0==strcmp(token,"vz"))
            vz = va_arg(args, double);
        if (0==strcmp(token,"a"))
            a = va_arg(args, double);
        if (0==strcmp(token,"P"))
            P = va_arg(args, double);
        if (0==strcmp(token,"e"))
            e = va_arg(args, double);
        if (0==strcmp(token,"inc"))
            inc = va_arg(args, double);
        if (0==strcmp(token,"Omega"))
            Omega = va_arg(args, double);
        if (0==strcmp(token,"omega"))
            omega = va_arg(args, double);
        if (0==strcmp(token,"pomega"))
            pomega = va_arg(args, double);
        if (0==strcmp(token,"f"))
            f = va_arg(args, double);
        if (0==strcmp(token,"M"))
            M = va_arg(args, double);
        if (0==strcmp(token,"E"))
            E = va_arg(args, double);
        if (0==strcmp(token,"l"))
            l = va_arg(args, double);
        if (0==strcmp(token,"theta"))
            theta = va_arg(args, double);
        if (0==strcmp(token,"T"))
            T = va_arg(args, double);
        if (0==strcmp(token,"h"))
            h = va_arg(args, double);
        if (0==strcmp(token,"k"))
            k = va_arg(args, double);
        if (0==strcmp(token,"ix"))
            ix = va_arg(args, double);
        if (0==strcmp(token,"iy"))
            iy = va_arg(args, double);
        if (0==strcmp(token,"primary")){
            primary = va_arg(args, struct reb_particle);
            primary_given = 1;
        }
        if (0==strcmp(token,"hash")){
            hash = va_arg(args, uint32_t);
        }
    }
    free(fmt_c);

    int Ncart = 0;
    if (!isnan(x)) Ncart++;
    if (!isnan(y)) Ncart++;
    if (!isnan(z)) Ncart++;
    if (!isnan(vx)) Ncart++;
    if (!isnan(vy)) Ncart++;
    if (!isnan(vz)) Ncart++;

    int Norb = 0;
    if (primary_given) Norb++;
    if (!isnan(a)) Norb++;
    if (!isnan(P)) Norb++;
    if (!isnan(e)) Norb++;
    if (!isnan(inc)) Norb++;
    if (!isnan(Omega)) Norb++;
    if (!isnan(omega)) Norb++;
    if (!isnan(pomega)) Norb++;
    if (!isnan(f)) Norb++;
    if (!isnan(M)) Norb++;
    if (!isnan(E)) Norb++;
    if (!isnan(l)) Norb++;
    if (!isnan(theta)) Norb++;
    if (!isnan(T)) Norb++;
    
    int Nnonpal = 0;
    if (primary_given) Nnonpal++;
    if (!isnan(e)) Nnonpal++;
    if (!isnan(inc)) Nnonpal++;
    if (!isnan(Omega)) Nnonpal++;
    if (!isnan(omega)) Nnonpal++;
    if (!isnan(pomega)) Nnonpal++;
    if (!isnan(f)) Nnonpal++;
    if (!isnan(M)) Nnonpal++;
    if (!isnan(E)) Nnonpal++;
    if (!isnan(theta)) Nnonpal++;
    if (!isnan(T)) Nnonpal++;
    
    int Npal = 0;
    if (!isnan(h)) Npal++;
    if (!isnan(k)) Npal++;
    if (!isnan(ix)) Npal++;
    if (!isnan(iy)) Npal++;
    
    int Nlong = 0;
    if (!isnan(f)) Nlong++;
    if (!isnan(M)) Nlong++;
    if (!isnan(E)) Nlong++;
    if (!isnan(l)) Nlong++;
    if (!isnan(theta)) Nlong++;
    if (!isnan(T)) Nlong++;

    if (Nnonpal>0 && Npal>0){
        *err = 7; // cannot mix pal and orbital elements
        return reb_particle_nan();
    }
    if (Ncart>0 && Norb>0){
        *err = 8; // cannot mix cartesian and orbital elements
        return reb_particle_nan();
    }

    if (Ncart || (!Norb)){ // Cartesian coordinates given, or not coordinates whatsoever
        struct reb_particle particle = {0};
        particle.hash = hash;
        particle.m = m;
        particle.r = radius;
        if (!isnan(x)) particle.x = x; // Note: is x is nan, then particle.x is 0  
        if (!isnan(y)) particle.y = y; 
        if (!isnan(z)) particle.z = z; 
        if (!isnan(vx)) particle.vx = vx; 
        if (!isnan(vy)) particle.vy = vy; 
        if (!isnan(vz)) particle.vz = vz; 
        return particle;
    }
    
    if (r==NULL){
        *err = 9; // need simulation for orbital elements
        return reb_particle_nan();
    }
    if (!primary_given){
#ifdef MPI
        reb_simulation_error(r, "When using MPI, you need to provide a primary to reb_simulation_add_fmt() when using orbital elements.");
        return reb_particle_nan();
#else // MPI
        primary = reb_simulation_com(r);
#endif // MPI
    }
    // Note: jacobi_masses not yet implemented.

    if (isnan(a) && isnan(P)){
        *err = 10; // can't have a and P
        return reb_particle_nan();
    }
    if (!isnan(a) && !isnan(P)){
        *err = 11; // need to have a or P
        return reb_particle_nan();
    }
    if (isnan(a)){
        a = cbrt(P*P*r->G *(primary.m + m)/(4.*M_PI*M_PI));
    }
    if (Npal>0){
        if (isnan(l)) l=0;
        if (isnan(h)) h=0;
        if (isnan(k)) k=0;
        if (isnan(ix)) ix=0;
        if (isnan(iy)) iy=0;
        if ((ix*ix + iy*iy) > 4.0){
            *err = 12; // e too high 
            return reb_particle_nan();
        }
        struct reb_particle particle = reb_particle_from_pal(r->G, primary, m, a, l, k, h, ix, iy);
        particle.r = radius;
        particle.hash = hash;
        return particle;
    }
    
    if (isnan(e)) e = 0.;
    if (isnan(inc)) inc = 0.;
    if (isnan(Omega)) Omega = 0.;
    
    if (!isnan(omega) && !isnan(pomega)){
        *err = 13; // Can't pass omega and pomega 
        return reb_particle_nan();
    }
    if (isnan(omega) && isnan(pomega)) omega = 0.;
    if (!isnan(pomega)){
        if (cos(inc)>0.){
            omega = pomega - Omega;
        }else{
            omega = Omega - pomega; // retrograde orbits
        }
    }

    if (Nlong>1){
        *err = 14; // only one longitude 
        return reb_particle_nan();
    }
    if (Nlong==0){
        f=0;
    }
    if (Nlong==1){
        if (!isnan(theta)){
            if (cos(inc)>0.){
                f = theta - Omega - omega;
            }else{
                f = Omega - omega - theta; // retrograde
            }
        }
        if (!isnan(l)){
            if (cos(inc)>0.){
                M = l - Omega - omega; // M will be converted to f below
            }else{
                M = Omega - omega - l; // retrograde
            }
        }
        if (!isnan(T)){
            double n = sqrt(r->G*(primary.m + m)/fabs(a*a*a));
            M = n * (r->t-T);
        }
        if (!isnan(M)){
            f = reb_M_to_f(e,M);
        }
        if (!isnan(E)){
            f = reb_E_to_f(e,E);
        }
    }
    struct reb_particle particle = reb_particle_from_orbit_err(r->G, primary, m, a, e, inc, Omega, omega, f, err);
    particle.r = radius;
    particle.hash = hash;
    return particle;
}

#define TINY 1.E-308 		///< Close to smallest representable floating point number, used for orbit calculation

struct reb_particle reb_particle_from_orbit_err(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f, int* err){
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

struct reb_particle reb_particle_from_orbit(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f){
	int err;
	return reb_particle_from_orbit_err(G, primary, m, a, e, inc, Omega, omega, f, &err);
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



// Calculates right quadrant for acos(num/denom) using a disambiguator that is < 0 when acos in the range (0, -pi)
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

struct reb_orbit reb_orbit_from_particle_err(double G, struct reb_particle p, struct reb_particle primary, int* err){
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

    hx = (dy*dvz - dz*dvy); 					// specific angular momentum vector
    hy = (dz*dvx - dx*dvz);
    hz = (dx*dvy - dy*dvx);
    o.h = sqrt ( hx*hx + hy*hy + hz*hz );		// abs value of angular momentum
    o.hvec.x = hx;
    o.hvec.y = hy;
    o.hvec.z = hz;

    vdiffsquared = vsquared - vcircsquared;	
    if(o.d <= TINY){							
        *err = 2;								// particle is on top of primary
        return reb_orbit_nan();
    }
    vr = (dx*dvx + dy*dvy + dz*dvz)/o.d;	
    rvr = o.d*vr;
    muinv = 1./mu;

    ex = muinv*( vdiffsquared*dx - rvr*dvx );
    ey = muinv*( vdiffsquared*dy - rvr*dvy );
    ez = muinv*( vdiffsquared*dz - rvr*dvz );
    o.e = sqrt( ex*ex + ey*ey + ez*ez );		// eccentricity
    o.evec.x = ex;
    o.evec.y = ey;
    o.evec.z = ez;
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

    double t0 = 0.0;
    if (p.sim != NULL){                     // if particle isn't in simulation yet, can't get time.
        t0 = p.sim->t;
    }
    o.T = t0 - o.M/fabs(o.n);               // time of pericenter passage (M = n(t-T).  Works for hyperbolic orbits using fabs and n as defined above).

    // move some of the angles into [0,2pi) range
    o.f = reb_mod2pi(o.f);
    o.l = reb_mod2pi(o.l);
    o.M = reb_mod2pi(o.M);
    o.theta = reb_mod2pi(o.theta);
    o.omega = reb_mod2pi(o.omega);
    
    
    // Cartesian eccentricity and inclination components, see Pal (2009)
    double fac = sqrt(2./(1.+hz/o.h))/o.h;
    o.pal_ix = -fac * hy;
    o.pal_iy = fac * hx;
    o.pal_k = o.h/mu*(dvy-dvz/(o.h+hz)*hy)-1./o.d*(dx-dz/(o.h+hz)*hx);
    o.pal_h = o.h/mu*(-dvx+dvz/(o.h+hz)*hx)-1./o.d*(dy-dz/(o.h+hz)*hy);
    return o;
}


struct reb_orbit reb_orbit_from_particle(double G, struct reb_particle p, struct reb_particle primary){
	int err;
	return reb_orbit_from_particle_err(G, p, primary, &err);
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
        double E = reb_M_to_E(e, M);
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

struct reb_particle reb_particle_from_pal(double G, struct reb_particle primary, double m, double a, double lambda, double k, double h, double ix, double iy){
    struct reb_particle np = {0};
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

void reb_simulation_rescale_var(struct reb_simulation* const r){
    // This function rescales variational particles if a coordinate
    // approached floating point limits (>1e100)
    if (r->N_var_config==0){
        return;
    }

    for (int v=0;v<r->N_var_config;v++){
        struct reb_variational_configuration* vc = &(r->var_config[v]);

        if (vc->lrescale <0 ) continue;  // Skip rescaling if lrescale set to -1

        int N = 1;
        if (vc->testparticle<0){
            N = r->N - r->N_var;
        }
        double scale = 0;
        struct reb_particle* const particles = r->particles + vc->index;
        for (int i=0; i<N; i++){
            struct reb_particle p = particles[i];
            scale = MAX(fabs(p.x), scale);
            scale = MAX(fabs(p.y), scale);
            scale = MAX(fabs(p.z), scale);
            scale = MAX(fabs(p.vx), scale);
            scale = MAX(fabs(p.vy), scale);
            scale = MAX(fabs(p.vz), scale);
        }
        if (scale > 1e100){
             
            if (vc->order == 1){
                for (int w=0;w<r->N_var_config;w++){
                    struct reb_variational_configuration* wc = &(r->var_config[w]);
                    if (wc->index_1st_order_a == vc->index || wc->index_1st_order_b == vc->index){
                        if (!(r->var_rescale_warning & 4)){
                            r->var_rescale_warning |= 4;
                            reb_simulation_warning(r, "Rescaling a set of variational equations of order 1 which are being used by a set of variational equations of order 2. Order 2 equations will no longer be valid.");
                        }
                    }
                }
            }else{ // order 2
                if (!(r->var_rescale_warning & 2)){
                    r->var_rescale_warning |= 2;
                    reb_simulation_warning(r, "Variational particles which are part of a second order variational equation have now large coordinates which might exceed range of floating point number range. REBOUND cannot rescale a second order variational equation as it is non-linear.");
                }
                return;
            }


            int is_synchronized = 1;
            if (r->integrator == REB_INTEGRATOR_WHFAST && r->ri_whfast.is_synchronized == 0){
                is_synchronized = 0;
            }
            if (r->integrator == REB_INTEGRATOR_EOS && r->ri_eos.is_synchronized == 0){
                is_synchronized = 0;
            }
            if (is_synchronized == 0){
                if (!(r->var_rescale_warning & 1)){
                    r->var_rescale_warning |= 1;
                    reb_simulation_warning(r, "Variational particles have large coordinates which might exceed range of floating point numbers. Rescaling failed because integrator was not synchronized. Turn on safe_mode or manually synchronize and rescale.");
                }
                return;
            }

            vc->lrescale += log(scale);
            for (int i=0; i<N; i++){
                particles[i].x /= scale;
                particles[i].y /= scale;
                particles[i].z /= scale;
                particles[i].vx /= scale;
                particles[i].vy /= scale;
                particles[i].vz /= scale;
            }

            if (r->integrator == REB_INTEGRATOR_WHFAST && r->ri_whfast.safe_mode == 0){
                r->ri_whfast.recalculate_coordinates_this_timestep = 1;
            }
        }
    }
}


int reb_simulation_add_variation_1st_order(struct reb_simulation* const r, int testparticle){
    r->N_var_config++;
    r->var_config = realloc(r->var_config,sizeof(struct reb_variational_configuration)*r->N_var_config);
    r->var_config[r->N_var_config-1].sim = r;
    r->var_config[r->N_var_config-1].order = 1;
    int index = r->N;
    r->var_config[r->N_var_config-1].index = index;
    r->var_config[r->N_var_config-1].lrescale = 0;
    r->var_config[r->N_var_config-1].testparticle = testparticle;
    struct reb_particle p0 = {0};
    if (testparticle>=0){
        reb_simulation_add(r,p0);
        r->N_var++;
    }else{
        int N_real = r->N - r->N_var;
        for (int i=0;i<N_real;i++){
            reb_simulation_add(r,p0);
        }
        r->N_var += N_real;
    }
    return index;
}


int reb_simulation_add_variation_2nd_order(struct reb_simulation* const r, int testparticle, int index_1st_order_a, int index_1st_order_b){
    r->N_var_config++;
    r->var_config = realloc(r->var_config,sizeof(struct reb_variational_configuration)*r->N_var_config);
    r->var_config[r->N_var_config-1].sim = r;
    r->var_config[r->N_var_config-1].order = 2;
    int index = r->N;
    r->var_config[r->N_var_config-1].index = index;
    r->var_config[r->N_var_config-1].lrescale = 0;
    r->var_config[r->N_var_config-1].testparticle = testparticle;
    r->var_config[r->N_var_config-1].index_1st_order_a = index_1st_order_a;
    r->var_config[r->N_var_config-1].index_1st_order_b = index_1st_order_b;
    struct reb_particle p0 = {0};
    if (testparticle>=0){
        reb_simulation_add(r,p0);
        r->N_var++;
    }else{
        int N_real = r->N - r->N_var;
        for (int i=0;i<N_real;i++){
            reb_simulation_add(r,p0);
        }
        r->N_var += N_real;
    }
    return index;
}

void reb_simulation_init_megno_seed(struct reb_simulation* const r, unsigned int seed){
	r->rand_seed = seed;
    reb_simulation_init_megno(r);
}

void reb_simulation_init_megno(struct reb_simulation* const r){
	r->megno_Ys = 0.;
	r->megno_Yss = 0.;
	r->megno_cov_Yt = 0.;
	r->megno_var_t = 0.;
	r->megno_n = 0;
	r->megno_mean_Y = 0;
	r->megno_mean_t = 0;
    int i = reb_simulation_add_variation_1st_order(r,-1);
	r->calculate_megno = i;
    const int imax = i + (r->N-r->N_var);
    struct reb_particle* const particles = r->particles;
    for (;i<imax;i++){ 
        particles[i].m  = 0.;
		particles[i].x  = reb_random_normal(r,1.);
		particles[i].y  = reb_random_normal(r,1.);
		particles[i].z  = reb_random_normal(r,1.);
		particles[i].vx = reb_random_normal(r,1.);
		particles[i].vy = reb_random_normal(r,1.);
		particles[i].vz = reb_random_normal(r,1.);
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
double reb_simulation_megno(struct reb_simulation* r){ // Returns the MEGNO <Y>
	if (r->t==0.) return 0.;
	return r->megno_Yss/r->t;
}
double reb_simulation_lyapunov(struct reb_simulation* r){ 
    // Returns the largest Lyapunov characteristic number (LCN)
    // Note that different definitions exist. 
    // Here, we're following Eq 24 of Cincotta and Simo (2000)
    // https://aas.aanda.org/articles/aas/abs/2000/20/h1686/h1686.html
	if (r->megno_var_t==0.0) return 0.;
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

void reb_tools_megno_update(struct reb_simulation* r, double dY, double dt_done){
	// Calculate running Y(t)
	r->megno_Ys += dY;
	double Y = r->megno_Ys/r->t;
	// Calculate averge <Y> 
	r->megno_Yss += Y * dt_done;
	// Update covariance of (Y,t) and variance of t
	r->megno_n++;
	double _d_t = r->t - r->megno_mean_t;
	r->megno_mean_t += _d_t/(double)r->megno_n;
	double _d_Y = reb_simulation_megno(r) - r->megno_mean_Y;
	r->megno_mean_Y += _d_Y/(double)r->megno_n;
	r->megno_cov_Yt += ((double)r->megno_n-1.)/(double)r->megno_n 
					*(r->t-r->megno_mean_t)
					*(reb_simulation_megno(r)-r->megno_mean_Y);
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

struct reb_vec3d reb_tools_spherical_to_xyz(const double magnitude, const double theta, const double phi){
    struct reb_vec3d xyz;
    xyz.x = magnitude * sin(theta) * cos(phi);
    xyz.y = magnitude * sin(theta) * sin(phi);
    xyz.z = magnitude * cos(theta);
    return xyz;
}  

void reb_tools_xyz_to_spherical(const struct reb_vec3d xyz, double* magnitude, double* theta, double* phi){
    *magnitude = sqrt(xyz.x*xyz.x + xyz.y*xyz.y + xyz.z*xyz.z);
    *theta = acos2(xyz.z, *magnitude, 1.);    // theta always in [0,pi] so pass dummy disambiguator=1
    *phi = atan2(xyz.y, xyz.x);
}  
