/**
 * @file 	librebound.c
 * @brief 	Declarations and functions for shared library.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 		Dave Spiegel <dave@ias.edu>
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Hanno Rein, Dave Spiegel, Daniel Tamayo
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
#include <string.h>
#include "particle.h"
#include "gravity.h"
#include "particle.h"
#include "main.h"
#include "librebound.h"
#include "boundaries.h"
#include "integrator.h"
#include "integrator_whfast.h"

// Default values of parameters and constants
double dt 	= 0.01;	
double t 	= 0;
double tmax	= 0;
double G 	= 1;
double softening = 0;
double timing 	= 0;
int closeEncounterPi = -1;
int closeEncounterPj = -1;
extern int Nmax;	
extern int iter;  // TODO DEBUG
const char *build_str = __DATE__ " " __TIME__;

// Function pointer to additional forces
void (*problem_additional_forces) (void) = NULL;
void (*problem_additional_forces_with_parameters) (struct particle* particles, double t, double dt, double G, int N, int N_megno) = NULL;
void (*problem_post_timestep_modifications) (void) = NULL;
void (*problem_post_timestep_modifications_with_parameters) (struct particle* particles, double t, double dt, double G, int N, int N_megno) = NULL;  

// Particle getter/setter methods.
void setp(struct particle* _p){
	free(particles);
	particles = malloc(sizeof(struct particle)*N);
	memcpy(particles,_p,sizeof(struct particle)*N);
}
struct particle particle_get(int i){
	return particles[i];
}
struct particle* particles_get(void){
	return particles;
}
void set_additional_forces(void (* _cb)(void)){
	problem_additional_forces = _cb;
}
void set_additional_forces_with_parameters(void (* _cb)(struct particle* particles, double t, double dt, double G, int N, int N_megno)){
	problem_additional_forces_with_parameters = _cb;
}
void set_post_timestep_modifications(void (* _cb)(void)){
	problem_post_timestep_modifications = _cb;
}
void set_post_timestep_modifications_with_parameters(void (* _cb)(struct particle* particles, double t, double dt, double G, int N, int N_megno)){
	problem_post_timestep_modifications_with_parameters = _cb;
}

void set_integrator(int i){
	integrator = i;
}

// Integrate for 1 step
void rebound_step(int do_timing){
    struct timeval tim;
	double timing_initial, timing_final;
	if (do_timing){
		gettimeofday(&tim, NULL);
		timing_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
	}
	integrator_part1();
	gravity_calculate_acceleration();
	if (N_megno){
		gravity_calculate_variational_acceleration();
	}
	if (problem_additional_forces) problem_additional_forces();
	if (problem_additional_forces_with_parameters) problem_additional_forces_with_parameters(particles, t, dt, G, N, N_megno);
	integrator_part2();
	if (problem_post_timestep_modifications){
		integrator_synchronize();
		problem_post_timestep_modifications();
		integrator_whfast_recalculate_jacobi_this_timestep = 1;
	}
	if (problem_post_timestep_modifications_with_parameters){
		integrator_synchronize();
		problem_post_timestep_modifications_with_parameters(particles, t, dt, G, N, N_megno);
		integrator_whfast_recalculate_jacobi_this_timestep = 1;
	}
	if (do_timing){
		gettimeofday(&tim, NULL);
		timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
		timing = timing_final-timing_initial;
	}
}

// Integrate until t=_tmax (or slightly more if exact_finish_time=0)
// Return values:
//   0 = All good
//   1 = No particles left
//   2 = Particle distance exceeds maxR
//   3 = Close encounter closer than minD
int rebound_integrate(double _tmax, int exact_finish_time, double maxR, double minD){
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double timing_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
	tmax = _tmax;
	double dt_last_done = dt;
	int last_step = 0;
	int ret_value = 0;
	const double dtsign = copysign(1.,dt); 				// Used to determine integration direction
	while(t*dtsign<tmax*dtsign && last_step<2 && ret_value==0){
		if (N<=0){
			fprintf(stderr,"\n\033[1mError!\033[0m No particles found. Exiting.\n");
			return(1);
		}
		rebound_step(0); 								// 0 to not do timing within step
		if ((t+dt)*dtsign>=tmax*dtsign && exact_finish_time==1){
			integrator_synchronize();
			dt = tmax-t;
			last_step++;
		}else{
			dt_last_done = dt;
		}
		if (maxR){
			// Check for escaping particles
			const double maxR2 = maxR*maxR;
			for (int i=0;i<N-N_megno;i++){
				struct particle p = particles[i];
				double r2 = p.x*p.x + p.y*p.y + p.z*p.z;
				if (r2>maxR2){
					ret_value = 2;
				}
			}
		}
		if (minD){
			// Check for close encounters
			const double minD2 = minD*minD;
			for (int i=0;i<N-N_megno;i++){
				struct particle pi = particles[i];
				for (int j=0;j<i;j++){
					struct particle pj = particles[j];
					const double x = pi.x-pj.x;
					const double y = pi.y-pj.y;
					const double z = pi.z-pj.z;
					const double r2 = x*x + y*y + z*z;
					if (r2<minD2){
						ret_value = 3;
						closeEncounterPi = i;
						closeEncounterPj = j;
					}
				}
			}
		}
	}
	integrator_synchronize();
	dt = dt_last_done;
	gettimeofday(&tim, NULL);
	double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	timing = timing_final-timing_initial;
	return ret_value;
}

void reset(void){ 
	dt 		= 0.01;
	t 		= 0;
	tmax		= 0;
	G 		= 1;
	softening 	= 0;
	N 		= 0;
	Nmax 		= 0;
	N_active 	= -1;
	N_megno 	= 0;
	iter		= 0;
	timing		= 0.;
	free(particles);
	particles 	= NULL;
	problem_additional_forces = NULL;
	problem_additional_forces_with_parameters = NULL;
	problem_post_timestep_modifications = NULL;
	problem_post_timestep_modifications_with_parameters = NULL;
	integrator_reset();
	struct timeval tim;
	gettimeofday(&tim, NULL);
	srand ( tim.tv_usec + getpid());
}
	 
